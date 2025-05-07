double theta_com_to_lab(double theta_com);

void draw_omp()
{
    bool setLogy = true;
    bool showRing = false;

    auto top_com = new LKDrawingGroup("top_com");
    auto top_lab = new LKDrawingGroup("top_lab");

    auto top1 = new LKDrawingGroup("top1");
    auto draw_a_com = top1 -> CreateDrawing("draw_a_com");
    auto draw_a_lab = top1 -> CreateDrawing("draw_a_lab");

    double legend_width = 0.45;
    draw_a_com -> SetCreateLegend(0,legend_width,0.07);
    draw_a_lab -> SetCreateLegend(1,legend_width,0.07);
    if (setLogy) {
        draw_a_com -> SetLogy();
        draw_a_lab -> SetLogy();
    }

    LKBinning bnnx_lab(100,0,70);
    LKBinning bnnx_com(100,25,180);
    LKBinning bnnx_com_x[3];
    bnnx_com_x[0] = LKBinning(100,40,80);
    bnnx_com_x[1] = LKBinning(100,100,120);
    bnnx_com_x[2] = LKBinning(100,126,162);
    LKBinning bnny(100,0,30);
    if (setLogy)
        bnny = LKBinning(100,0.1,500);

    auto hist_com = (bnnx_com*bnny).NewH2("hist_com",";#theta_{com};Cross section");
    hist_com -> SetStats(0);
    draw_a_com -> Add(hist_com,"",".");

    auto hist_lab = (bnnx_lab*bnny).NewH2("hist_lab",";#theta_{lab};Cross section");
    hist_lab -> SetStats(0);
    draw_a_lab -> Add(hist_lab,"",".");

    for (auto hist : {hist_com,hist_lab}) {
        hist -> GetXaxis() -> SetLabelFont(132);
        hist -> GetYaxis() -> SetLabelFont(132);
        hist -> GetXaxis() -> SetTitleFont(132);
        hist -> GetYaxis() -> SetTitleFont(132);
        hist -> GetXaxis() -> SetLabelSize(0.045);
        hist -> GetYaxis() -> SetLabelSize(0.045);
        hist -> GetXaxis() -> SetTitleSize(0.050);
        hist -> GetYaxis() -> SetTitleSize(0.050);
        hist -> GetXaxis() -> SetTitleOffset(1.12);
        hist -> GetYaxis() -> SetTitleOffset(1.2);
    }

    map<int, double> df;
    df[21] = 0.505;
    df[23] = 0.487;
    df[24] = 0.377;
    df[25] = 0.310;
    //vector<int> NaArray = {21,23,24,25};
    vector<int> NaArray = {21,25};

    LKParameter colors("colors.mac");
    if (NaArray.size()==2) colors = LKParameter("colors2.mac");

    int count = 0;
    for (auto ANa : NaArray)
    {
        auto draw_com = top_com -> CreateDrawing("com",true);
        auto draw_lab = top_lab -> CreateDrawing("lab",true);
        draw_lab -> SetCanvasMargin(0.12,0.05,0.13,0.05);
        draw_com -> SetCanvasMargin(0.12,0.05,0.13,0.05);
        draw_lab -> Add(hist_lab,".");
        draw_com -> Add(hist_com,".");
        draw_com -> AddLegendLine(Form("#beta_{2} = %.3f",df[ANa]));
        //draw_lab -> AddLegendLine(Form("#beta_{2} = %.3f",df[ANa]));
        draw_com -> SetCreateLegend(0,0.40,0.08);
        draw_lab -> SetCreateLegend(1,0.40,0.08);
        if (setLogy) {
            draw_com -> SetLogy();
            draw_lab -> SetLogy();
        }

        for (auto iType : {0,1})
        {
            TString fileNameOMP = Form("data_omp/cx_avg_drhbc_%dNa.txt",ANa);
            if (iType==1) fileNameOMP = Form("data_omp/cx_rchb_%dNa.txt",ANa);
            auto graph_com = new TGraph();
            auto graph_lab = new TGraph();
            ifstream file(fileNameOMP);
            std::string line;
            std::getline(file, line);
            std::getline(file, line);
            double theta_com, y;

            while (file >> theta_com >> y)
            {
                if ((bnnx_com*bnny).IsInside(theta_com,y))
                {
                    auto theta_lab = theta_com_to_lab(theta_com);
                    graph_com -> SetPoint(graph_com->GetN(),theta_com,y);
                    graph_lab -> SetPoint(graph_lab->GetN(),theta_lab,y);
                    /*
                    if (iType==0)
                    {
                        for (auto bnnx : {bnnx_com1,bnnx_com2,bnnx_com3})
                        {
                            if (bnnx.IsInside(theta_com))
                            {
                                double ee = 0.1;
                                double yh = y + y * (TMath::Power(10,ee) - 1);
                                double yl = y - y * (1 - TMath::Power(10,-ee));
                                graph_com_s -> SetPoint(graph_com_s->GetN(),theta_com,y);
                                graph_com_s -> SetPointError(graph_com_s->GetN()-1, 0, 0, yl, yh);
                                cout << ANa << " " << theta_com << " " << yl << " " << yh << endl;
                            }
                            draw_com -> Add(graph_com_s,"3",".");
                        }
                    }
                    */
                }
            }
            file.close();

            auto title_reaction = Form("^{%d}Na+p",ANa);
            //auto title = Form("Deformed ^{%d}Na+p",ANa);
            //if (iType==1) title = Form("Spherical ^{%d}Na+p",ANa);
            auto title = "Deformed";
            if (iType==1) title = "Spherical";
            graph_com -> SetLineColor(colors.GetColor(count));
            graph_lab -> SetLineColor(colors.GetColor(count));
            graph_com -> SetLineWidth(2);
            graph_lab -> SetLineWidth(2);
            if (iType==1) {
                graph_com -> SetLineStyle(2);
                graph_lab -> SetLineStyle(2);
            }

            draw_a_com -> Add(graph_com,"l",title);
            draw_a_lab -> Add(graph_lab,"l",title);
            draw_com -> Add(graph_com,"l",title);
            draw_lab -> Add(graph_lab,"l",title);

            //auto tt = new TLatex(bnnx_com.Lerp(1.-legend_width),bnny.Lerp(0.95),title_reaction);
            //tt -> SetTextAlign(33);
            //tt -> SetTextFont(132);
            //tt -> SetTextSize(0.08);
            //draw_com -> Add(tt);

            //auto pv = new TPaveText(bnnx_com.x1()+0.03*bnnx_com.d(),bnny.x2()-0.18*bnny.d(),bnnx_com.x1()+0.3*bnnx_com.d(),bnny.x2());
            //auto pv = new TPaveText(bnnx_com.x1()+0.03*bnnx_com.d(),bnny.x2()-0.18*bnny.d(),bnnx_com.x1()+0.3*bnnx_com.d(),bnny.x2());
            //auto pv = new TPaveText(0.03,0.82,0.3,1,"NDC");
            auto pv = new TPaveText(0.15,0.80,0.4,0.92,"NDC");
            pv -> AddText(title_reaction);
            pv -> SetTextAlign(22);
            pv -> SetTextFont(132);
            pv -> SetTextSize(0.08);
            pv -> SetBorderSize(1);
            draw_com -> Add(pv,"same auto");
        }

        if (0)
        if (ANa==23)
        {
            TString fileNamePre = Form("data_pre/x_%dNapp_8MeV.csv",ANa);
            auto graph = new TGraph(fileNamePre,"%lg, %lg");
            graph -> SetMarkerStyle(20);
            graph -> Print();
            draw_com -> Add(graph,"p","J. Hellstrom (1969)");
        }

        count++;
    }

    //top1 -> SetCanvasSize(1200,500);
    //top1 -> Draw();

    for (auto top : {top_com,top_lab}) { 
        //top -> SetCanvasDivision(4,1);
        top -> SetCanvasDivision(2,2);
        if (NaArray.size()==2)
            top -> SetCanvasDivision(1,2);
        top -> Draw();
        top -> GetCanvas() -> SaveAs(Form("figures/fig_%s_%d%d%s.png",top->GetName(),top->GetDivX(),top->GetDivY(),(setLogy?".logy":""),(showRing?".ring":"")));
    }
}

double theta_com_to_lab(double theta_com) { return (180. - theta_com)/2.; }
