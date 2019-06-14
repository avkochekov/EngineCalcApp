#include "mainwindow.h"
#include "kinematic.h"
#include "ui_mainwindow.h"
#include "math.h"
#include <QDialog>

    double d;
    double d_g;
    double l_sh;
    double h_g;
    double s_v;

    double d_2sh;
    double C_b;
    double t_v;
    double l_k;

    double h_sh_min;
    double h_sh;
    double b_sh;
    double a_sh;
    double t_sh;

    double n_N;
    double E_sh;
    double E_v;
    double Eps_m;
    double Eps_p;


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->centralWidget->setMaximumSize(760, 560);
    ui->centralWidget->setMinimumSize(760, 560);

    ui->tabWidget_4->removeTab(3);
//    ui->tabWidget_4->removeTab(4);

    // Вывод варнинга о нефункциональности проги...
    QDate   currentDate = QDate::currentDate();
    QDate   limitDate(2015, 9, 01);

    if (currentDate >= limitDate){
        QString InfoTitle = "Информация";
        QString InfoText = "Функционал данной програмы по определенным причинам реализован не полностью. "
                "В рассчетах возможны неточночти и ошибки."
                "\n\n"
                "По вопросам доработки программы можете обратиться на \ne-mail: crossprog55@gmail.com";
        QMessageBox InfoMessage;
        InfoMessage.information(ui->centralWidget, InfoTitle, InfoText, "Ok");
    }

    QString AboutProgramm = "Программа расчета параметров двигателя.\n"
            "CrossProg55@gmail.com\n"
            "Qt " + QString(QT_VERSION_STR);
    ui->lbl_AboutProgramm->setText(AboutProgramm);
    // При запуске программы активной будет первая вкладка (ui->tab)
    ui->tabWidget->setCurrentWidget(ui->tab);
}

MainWindow::~MainWindow()
{
    delete ui;
}

//void MainWindow::SetTableVision(QTableWidget table, int ccount, int rcount, QStringList hlabels)
//{
//    table->setColumnCount(ccount);
//    table->setRowCount(rcount);
//    table->setHorizontalHeaderLabels(hlabels);
//}

//void MainWindow::PlotGraphic(QCustomPlot plot, QVector xData, QVector yData, QString xLabel, QString yLabel, int yMin, int yMax){
//    plot->clearGraphs();
//    plot->addGraph();
//    plot->graph()->setData(xData,yData);
//    plot->graph()->setLineStyle((QCPGraph::lsLine));
//    plot->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
//    plot->xAxis->setLabel(xLabel);
//    plot->yAxis->setLabel(yLabel);
//    plot->xAxis->setRange(0, 720/2);
//    plot->yAxis->setRange(1.05*yMin, 1.05*yMax);
//    plot->replot();
//}

void MainWindow::on__btn_CalcKinematic_clicked()
{
    double R = ui->param_R->text().toDouble();
    double L = ui->param_L->text().toDouble();
    double lambda = R/L;
        ui->param_lambda->setText(QString::number(lambda));

    double n_N = ui->param_n_N->text().toDouble();
    double omega = 3.14 * n_N / 30;
        ui->param_omega->setText(QString::number(omega));

    QVector<double> Sx(73), C_phi(73), J_phi(73);
    QVector<double> Deg(73);

    QStringList HeaderLabels;
    HeaderLabels << "φ(deg)"
                 << "S(φ)"
                 << "C(φ)"
                 << "J(φ)";

//    SetTableVision(ui->table, 4,73,HeaderLabels);
    ui->table->setColumnCount(4);
    ui->table->setRowCount(73);
    ui->table->setHorizontalHeaderLabels(HeaderLabels);

    double maxY_Sx = 0, maxY_C_phi = 0, maxY_J_phi = 0;
    double minY_Sx = 0, minY_C_phi = 0, minY_J_phi = 0;

    for (int row=0; row<ui->table->rowCount(); row++){
        double deg = M_PI / 180;
        double phi = row * 10 * deg;

        Deg[row] = row * 10;
        Sx[row]      = R * ((1-cos(phi)) + 1/4 * (1-cos(2 * phi)));
            if (Sx[row] > maxY_Sx) maxY_Sx = Sx[row];
            if (Sx[row] < minY_Sx) minY_Sx = Sx[row];
        C_phi[row]   = R * omega * (sin(phi) + lambda/2 * sin(2*phi));
            if (C_phi[row] > maxY_C_phi) maxY_C_phi = C_phi[row];
            if (C_phi[row] < minY_C_phi) minY_C_phi = C_phi[row];
        J_phi[row]   = R * pow(omega,2) * (cos(phi) + lambda * cos(2*phi));
            if (J_phi[row] > maxY_J_phi) maxY_J_phi = J_phi[row];
            if (J_phi[row] < minY_J_phi) minY_J_phi = J_phi[row];

        for (int col=0; col<ui->table->columnCount(); col++){
            QTableWidgetItem *Item = new QTableWidgetItem();
            double VarItem;
            switch (col) {
            case 0:
                VarItem = row*10;
                break;
            case 1:
                VarItem = Sx[row];
                break;
            case 2:
                VarItem = C_phi[row];
                break;
            case 3:
                VarItem = J_phi[row];
                break;
            default:
                break;
            }
            if (col == 0)
                Item->setText(QString::number((int)VarItem));
            else
                Item->setText(QString::number(VarItem, 'f', 4));
            Item->setTextAlignment(Qt::AlignCenter);
            ui->table->setItem(row,col,Item);
            ui->table->setColumnWidth(col,100);
        }
    }

//    PlotGraphic(ui->plot_porsh_position,    Deg, Sx,    "Degree", "Sx",     minY_Sx, maxY_Sx);
//    PlotGraphic(ui->plot_porsh_speed,       Deg, C_phi, "Degree", "C(phi)", minY_C_phi, maxY_C_phi);
//    PlotGraphic(ui->plot_porsh_accel,       Deg, J_phi, "Degree", "J(phi)", minY_J_phi, maxY_J_phi);

    ui->plot_porsh_position->clearGraphs();
    ui->plot_porsh_position->addGraph();
    ui->plot_porsh_position->graph()->setData(Deg,Sx);
    ui->plot_porsh_position->graph()->setLineStyle((QCPGraph::lsLine));
    ui->plot_porsh_position->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
//    ui->plot_widget->plotLayout()->insertRow(0);
//    ui->plot_widget->plotLayout()->addElement(0, 0, new QCPPlotTitle(ui->plot_widget, "Абра-кадабра"));
    ui->plot_porsh_position->xAxis->setLabel("φ");
    ui->plot_porsh_position->yAxis->setLabel("Sx(φ)");
    ui->plot_porsh_position->xAxis->setRange(0, 720/2);
    ui->plot_porsh_position->yAxis->setRange(1.05*minY_Sx, 1.05*maxY_Sx);
    ui->plot_porsh_position->replot();

    ui->plot_porsh_speed->clearGraphs();
    ui->plot_porsh_speed->addGraph();
    ui->plot_porsh_speed->graph()->setData(Deg,C_phi);
    ui->plot_porsh_speed->graph()->setLineStyle((QCPGraph::lsLine));
    ui->plot_porsh_speed->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
    ui->plot_porsh_speed->xAxis->setLabel("φ");
    ui->plot_porsh_speed->yAxis->setLabel("C(φ)");
    ui->plot_porsh_speed->xAxis->setRange(0, 720/2);
    ui->plot_porsh_speed->yAxis->setRange(1.05*minY_C_phi, 1.05*maxY_C_phi);
    ui->plot_porsh_speed->replot();

    ui->plot_porsh_accel->clearGraphs();
    ui->plot_porsh_accel->addGraph();
    ui->plot_porsh_accel->graph()->setData(Deg,J_phi);
    ui->plot_porsh_accel->graph()->setLineStyle((QCPGraph::lsLine));
    ui->plot_porsh_accel->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
    ui->plot_porsh_accel->xAxis->setLabel("φ");
    ui->plot_porsh_accel->yAxis->setLabel("J(φ)");
    ui->plot_porsh_accel->xAxis->setRange(0,720/2);
    ui->plot_porsh_accel->yAxis->setRange(1.05*minY_J_phi, 1.05*maxY_J_phi);
    ui->plot_porsh_accel->replot();
}
void MainWindow::on__btn_CalcTemp_clicked()
{
    double C =  ui->param_temp_C->text().toDouble();
    double H =  ui->param_temp_H->text().toDouble();
    double O =  ui->param_temp_O->text().toDouble();
    double a =  ui->param_temp_a->text().toInt();
    double e =  ui->param_temp_e->text().toDouble();
    double ez = ui->param_temp_E->text().toDouble();
    double vbp =ui->param_temp_vbp->text().toInt();
    double la = ui->param_temp_lambda->text().toDouble();
    double t1 = ui->param_temp_t1->text().toDouble();
    double a1 = ui->param_temp_A->text().toDouble();
    double b =  ui->param_temp_B->text().toDouble();
    double ne = ui->param_temp_Ne->text().toDouble();
    double n =  ui->param_temp_n->text().toDouble();

    double l, l1, m1;
    double mh2o, mco2, mn2;
    double mh2, m2, pk1, po, pa, pa1;
    double t, tr, pr, mr, yr, ta, Nv;
    double n1, pc, tc, tc1, mcv, mcvco2, mcvco, mcvh2o, mcvh2, mcvn2, mcvo2;
    double mcv2, mcv1, my0, my, hu1, hu, hrab, mcv2t;
    double mcp, tz, pz, ro, tz1;
    double n2, k2, si, tb, tr1, pi1, PI, fi, nyi, gi, pm, Pe, nym, nye, ge, gt;
    double vncr, vl, Vh, d, S, fn, me;
    double pb;

    double K = 0.461;
    double pk = 0.17;
    double tk = 293;
    double rb = 287;

    double mo2, mco;
    double mcp2;



    l = (1 / 0.21) * ((C / 12) + (H / 4) - (O / 32));
        l1 = (1 / 0.23) * ((8 / 3 * C) + (8 * H) - O);
        m1 = a * l;

        if (a >= 1){
            mh2o = H / 2;
            mco2 = C / 12;
            mn2 = (0.79 * a) * l;
            mo2 = 0.21 * (a - 1) * l;
        } else {
            mco2 = (C / 12) - 2 * ((1 - a) / (1 + K)) * 0.21 * l;
            mco = 2 * ((1 - a) / (1 + K)) * 0.21 * l;
            mh2o = (H / 2) - 2 * K * ((1 - a) / (1 + K)) * 0.21 * l;
            mh2 = 2 * K * ((1 - a) / (1 + K)) * 0.21 * l;
            mn2 = 0.79 * a * l;
        }

        if (a > 1){
            m2 = mco2 + mh2o + mn2 + mo2;
        } else if (a = 1){
            m2 = mco2 + mh2o + mn2;
        } else {
            m2 = mco2 + mco + mh2o + mh2 + mn2;
        }

        pk1 = (pk * (pow(10,6) / (rb * tk)));
        po = 3 * (pow((vbp),2) / 2) * pk * 0.000001;
        pa = pk - pa1;
        pr = 0.115;
        t = 25;
        tr = 700;
        tk = 293;
        yr = (((tk + t) / tr) * ((pr) / (e * po - pr)));
        mr = yr * m1;
        ta = (tk + t + yr * tr) / (1 + yr);
        Nv = ((tk) / (tk + t)) * ((1) / (e - 1)) * (1 / pk) * (e * pa - pr);
        n1 = 1.35;
        pc = pa * pow(e,n1);
        tc = ta * pow(e,(n1 - 1));
        tc1 = tc - 273;
        mcv = 20.6 + 2.638 * 0.001 * tc1;
        mcvco2 = 27.941 + 0.019 * t1 - 0.000005487 * pow(t1,2);
        mcvco = 20.597 + 0.00267 * t1;
        mcvh2o = 24.953 + 0.005359 * t1;
        mcvh2 = 20.684 + 0.000206 * t1 + 0.000000588 * pow(t1,2);
        mcvn2 = 20.398 + 0.0025 * t1;
        mcvo2 = 20.93 + 0.0046441 * t1 - 0.00000084 * pow(t1,2);
        mcv2 = (1 / m2) * (mco2 * mcvco2 + mco * mcvco + mh2o * mcvh2o + mh2 * mcvh2 + mn2 * mcvn2 + mo2 * mcvo2);
        mcv1 = (1 / (1 + yr)) * (mcv + yr * mcv2);
        my0 = m2 / m1;
        my = (my0 + yr) / (1 + yr);

        if (a < 1) hu1 = 119950 * (1 - a) * l;
        hu = 42440;
        hrab = (hu - hu1) / (m1 * (1 + yr));

        if (a >= 1)
            mcv2t = (1 / m2) * (mco2 * (mcvco2) + mh2o * (mcvh2o) + mn2 * (mcvn2) + mo2 * (mcvo2));
        else
            mcv2t = (1 / m2) * (mco2 * (mcvco2) + mco * (mcvco) + mh2o * (mcvh2o) + mh2 * (mcvh2) + mn2 * (mcvn2));

        mcp2 = mcv2t + 8.315;
        tz = (ez * hrab + (mcv1 + 8.315 * la) * tc1 + 2270 * (la - my)) / (my * mcp2);
        pz = la * pc;
        tz1 = tz + 273;
        ro = (my * tz1) / (la * tc);
        n2 = k2 = 1.25;
        si = e / ro;
        pb = pz / pow(si,n2);
        tb = tz / pow(si,(n2 - 1));
        tr1 = tb / pow((sqrt(pb / pr)),(1 / 3));
        pi1 = (pc / (e - 1)) * (((la * ro) / (n2 - 1)) * (1 - (1 / pow(si,(n2 - 1)))) + la * (ro - 1) - (1 / (n1 - 1)) * (1 - (1 / pow(e,(n1 - 1)))));
        fi = 0.93;
        PI = fi * pi1;

        nyi = (PI * l * a) / (42440 * ro * Nv);

        gi = (3600 / (42440 * nyi));
        pm = a1 + b * vncr;
        Pe = PI - pm;
        nym = Pe / PI;
        nye = nym * nyi;
        ge = 3600 / (hu * nye);
        gt = ge * ne * 0.001;
            // Изначально vl = (30 * 4 * ne) / (Pe * n);
        vl = (30 * 4 * ne) / (Pe * n);
        Vh = vl / 8;
        d = 100 * pow((sqrt((4 * Vh) / (3.14 * 1))),(1 / 3));
        fn = (3.14 * pow(d,2)) / 4;
            // Изначально me = (30000 / 3.14) * (ne / n);
        me = (30000 / 3.14) * (ne / n);
        vncr = (d * ne) / 30000;

        ui->param_temp_L0->setText(QString::number(l, 'f', 4));
        ui->param_temp_l0->setText(QString::number(l1, 'f', 4));
        ui->param_temp_M1->setText(QString::number(m1, 'f', 4));
        ui->param_temp_M2->setText(QString::number(m2, 'f', 4));
        ui->param_temp_pk->setText(QString::number(pk, 'f', 4));
        ui->param_temp_P0->setText(QString::number(po, 'f', 4));
        ui->param_temp_Pa->setText(QString::number(pa, 'f', 4));
        ui->param_temp_y->setText(QString::number(yr, 'f', 4));
        ui->param_temp_Mr->setText(QString::number(mr, 'f', 4));
        ui->param_temp_Ta->setText(QString::number(ta, 'f', 4));
        ui->param_temp_nv->setText(QString::number(Nv, 'f', 4));
        ui->param_temp_Pc->setText(QString::number(pc, 'f', 4));
        ui->param_temp_Tc->setText(QString::number(tc, 'f', 4));
        ui->param_temp_mcv->setText(QString::number(mcv, 'f', 4));
        ui->param_temp_mcv2->setText(QString::number(mcv2, 'f', 4));
        ui->param_temp_mcv1->setText(QString::number(mcv1, 'f', 4));
        ui->param_temp_my0->setText(QString::number(my0, 'f', 4));
        ui->param_temp_my->setText(QString::number(my, 'f', 4));
        ui->param_temp_hu1->setText(QString::number(hu1, 'f', 4));
        ui->param_temp_hrab->setText(QString::number(hrab, 'f', 4));
        ui->param_temp_mcv2t->setText(QString::number(mcv2t, 'f', 4));
        ui->param_temp_mcp->setText(QString::number(mcp, 'f', 4));
        ui->param_temp_tz->setText(QString::number(tz, 'f', 4));
        ui->param_temp_pz->setText(QString::number(pz, 'f', 4));
        ui->param_temp_ro->setText(QString::number(ro, 'f', 4));
        ui->param_temp_tz1->setText(QString::number(tz1, 'f', 4));
        ui->param_temp_si->setText(QString::number(si, 'f', 4));
        ui->param_temp_pb->setText(QString::number(pb, 'f', 4));
        ui->param_temp_tb->setText(QString::number(tb, 'f', 4));
        ui->param_temp_tr->setText(QString::number(tr, 'f', 4));
        ui->param_temp_pi1->setText(QString::number(pi1, 'f', 4));
        ui->param_temp_pi->setText(QString::number(PI, 'f', 4));
        ui->param_temp_nyi->setText(QString::number(nyi, 'f', 4));
        ui->param_temp_gi->setText(QString::number(gi, 'f', 4));
        ui->param_temp_pm->setText(QString::number(pm, 'f', 4));
        ui->param_temp_pe->setText(QString::number(Pe, 'f', 4));
        ui->param_temp_fn->setText(QString::number(fn, 'f', 4));
        ui->param_temp_nym->setText(QString::number(nym, 'f', 4));
        ui->param_temp_nye->setText(QString::number(nye, 'f', 4));
        ui->param_temp_ge->setText(QString::number(ge, 'f', 4));
        ui->param_temp_gt->setText(QString::number(gt, 'f', 4));
        ui->param_temp_vl->setText(QString::number(vl, 'f', 4));
        ui->param_temp_vh->setText(QString::number(Vh, 'f', 4));
        ui->param_temp_d->setText(QString::number(d, 'f', 4));
        ui->param_temp_S->setText(QString::number(S, 'f', 4));
        ui->param_temp_me->setText(QString::number(me, 'f', 4));
        ui->param_temp_vncr->setText(QString::number(vncr, 'f', 4));
}
void MainWindow::on__btn_CalcEng_clicked()
{
    double Pe   = ui->param_eng_Pe      ->text().toDouble();
    double i    = ui->param_eng_i       ->text().toDouble();
    double Vh   = ui->param_eng_Vh      ->text().toDouble();
    double Vl   = ui->param_eng_Vl      ->text().toDouble();
    double nNorm= ui->param_eng_nNorm   ->text().toDouble();
    double nMin = ui->param_eng_nMin    ->text().toDouble();
    double Dp   = ui->param_eng_Dp      ->text().toDouble();
    double S    = ui->param_eng_S       ->text().toDouble();
    double M    = ui->param_eng_M       ->text().toDouble();
    double tau  = ui->param_eng_tau     ->text().toDouble();

    // Эффективная мощность
    double Nv;
    double NvMin;
        Nv = (Pe * Vh * i * nNorm) / (30 * tau);
        NvMin = (Pe * Vh * i * nMin) / (30 * tau);
    ui->param_eng_Nv->setText(QString::number(Nv, 'f', 4));
    ui->param_eng_NvMin->setText(QString::number(NvMin, 'f', 4));

    // Диаметр цилиндра 1
    double Dc;
    double DcMin;
        Dc = pow(((120 * tau * Nv) / (M_PI * Pe * nNorm * i * (S / Dp))),(1 / 3));
        DcMin = pow(((120 * tau * NvMin) / (M_PI * Pe * nMin * i * (S / Dp))),(1 / 3));
    ui->param_eng_Dc->setText(QString::number(Dc, 'f', 4));
    ui->param_eng_DcMin->setText(QString::number(DcMin, 'f', 4));

    // Литрвая мощность
    double Nl;
    double NlMin;
        Nl = (Pe * nNorm) / (30 * tau);
        NlMin = (Pe * nMin) / (30 * tau);
    ui->param_eng_Nl->setText(QString::number(Nl, 'f', 4));
    ui->param_eng_NlMin->setText(QString::number(NlMin, 'f', 4));

    // Диаметр цилиндра 2
    double DcWithNl;
    double DcWithNlMin;
    DcWithNl = pow(((4 * Nv) / (M_PI * Nl * i * (S / Dp))),(1 / 3));
    DcWithNlMin = pow(((4 * Nv) / (M_PI * NlMin * i * (S / Dp))),(1 / 3));
    ui->param_eng_DcNl->setText(QString::number(DcWithNl, 'f', 4));
    ui->param_eng_DcNlMin->setText(QString::number(DcWithNlMin, 'f', 4));

    // Оценка быстроходности
    double Cm;
    double CmMin;
    Cm = (S * nNorm) / 30;
    CmMin = (S * nMin) / 30;
    ui->param_eng_Cm->setText(QString::number(Cm, 'f', 4));
    ui->param_eng_CmMin->setText(QString::number(CmMin, 'f', 4));

    // Удельная масса двигателя
    double mN;
    double mNMin;
    mN = M / Nv;
    mNMin = M / NvMin;
    ui->param_eng_mN->setText(QString::number(mN, 'f', 4));
    ui->param_eng_mNMin->setText(QString::number(mNMin, 'f', 4));

    // Отношение хода поршня к диаметру цилиндра
    ui->param_eng_SD->setText(QString::number(S/Dc, 'f', 4));
    ui->param_eng_SDMin->setText(QString::number(S/DcMin, 'f', 4));
}
void MainWindow::on__btn_CalcGroup0_Apply_clicked()
{
    d       = ui->param_Group0_d    ->text().toDouble();
    d_g     = ui->param_Group0_dg   ->text().toDouble();
    l_sh    = ui->param_Group0_lsh  ->text().toDouble();
    h_g     = ui->param_Group0_hg   ->text().toDouble();
    s_v     = ui->param_Group0_sv   ->text().toDouble();

    d_2sh   = ui->param_Group0_d2sh ->text().toDouble();
    C_b     = ui->param_Group0_cb   ->text().toDouble();
    t_v     = ui->param_Group0_tv   ->text().toDouble();
    l_k     = ui->param_Group0_lk   ->text().toDouble();

    h_sh_min= ui->param_Group0_hshMin->text().toDouble();
    h_sh    = ui->param_Group0_hsh  ->text().toDouble();
    b_sh    = ui->param_Group0_bsh  ->text().toDouble();
    a_sh    = ui->param_Group0_ash  ->text().toDouble();
    t_sh    = a_sh;

    n_N     = ui->param_Group0_nN   ->text().toDouble();
    E_sh    = ui->param_Group0_Esh  ->text().toDouble();
    E_v     = ui->param_Group0_Ev   ->text().toDouble();
    Eps_m   = ui->param_Group0_Epsm ->text().toDouble();
    Eps_p   = ui->param_Group0_Epsp ->text().toDouble();

}
void MainWindow::on__btn_CalcGroup_1_clicked()
{
    double m_sh;
    double Sigma_v;
    double Delta;
    double dT;
    double d_n;

    double Pj, m_vg, Omg_xx_max, Sigma_max, k_sigma, Sigma_mo, Sigma_ak0, Delta_t, Delta_E, p, Sigma_ap_d, Sigma_ap_i;

    // Неизвестные переменные
    double n_xx_max, omg_v, m_p, R, phi, lambda, Sigma_a0, eps_M, eps_N, alpha_v, alpha_m, p1, p2;

    m_sh    = ui->param_Group1_msh    ->text().toDouble();
    Sigma_v = ui->param_Group1_Sigmav ->text().toDouble();
    Delta   = ui->param_Group1_Delta  ->text().toDouble();
    dT      = ui->param_Group1_DeltaT ->text().toDouble();
    d_n     = ui->param_Group1_dn     ->text().toDouble();

    //  Поршневая головка шатуна
    //   Расчет сечения I - I

    // Максимальное напряжение пульсирующего цикла
    m_vg = 0.06 * m_sh;                              // m_sh      ?
    Omg_xx_max = 3.14 * n_xx_max / 30;               // n_xx_max  ?
    k_sigma = 1.2 + 1.8 * 0.0001 * (omg_v - 400);    // omg_v     ?

    // Переменная сила инерции
    Pj = -(m_p + m_vg) * pow(Omg_xx_max,2) * R * (cos(phi) + lambda * cos(2 * phi));

    // Выбор eps_m и  eps_n

        Sigma_max = ((m_p + m_vg) * pow(Omg_xx_max,2) * R * (1 + lambda)) / 2 * h_g * l_sh;     // m_p, R, lambda    ?
        Sigma_mo = Sigma_max / 2;
        Sigma_a0 = Sigma_mo;
        Sigma_ak0 = Sigma_a0 * k_sigma / (eps_M * eps_N);

    // Расчет напряжения от запрессованой втулки
        Delta_t = d * (alpha_v + alpha_m) * dT;      // alpha_v, alpha_m  ?
        Delta_E = Delta + Delta_t;                   // Delta             ?

        p1 = (pow(d_g,2) + pow(d,2)) / (pow(d_g,2) - pow(d,2)) + 0.3;
        p2 = (pow(d,2) + pow(d_n,2)) / (pow(d,2) - pow(d_n,2)) - 0.3;
        p = Delta_E / (d * (p1 / E_sh + p2 / E_v));  // E_sh, E_v         ?

        Sigma_ap_d = p * 2 * pow(d,2) / (pow(d_g,2) - pow(d,2));
        Sigma_ap_i = p * (pow(d_g,2) + pow(d,2)) / (pow(d_g,2) - pow(d,2));

//    ui->param_Group1_
    ui->param_Group1_Pjp        ->setText(QString::number(Pj, 'f', 4));
    ui->param_Group1_mvg        ->setText(QString::number(m_vg, 'f', 4));
    ui->param_Group1_omg_xx_max ->setText(QString::number(Omg_xx_max, 'f', 4));
    ui->param_Group1_SigmaMax   ->setText(QString::number(Sigma_max, 'f', 4));
    ui->param_Group1_ksigma     ->setText(QString::number(k_sigma, 'f', 4));
    ui->param_Group1_Sigmam0    ->setText(QString::number(Sigma_mo, 'f', 4));
    ui->param_Group1_Sigmaak    ->setText(QString::number(Sigma_ak0, 'f', 4));
    ui->param_Group1_Deltat     ->setText(QString::number(Delta_t, 'f', 4));
    ui->param_Group1_DeltaE     ->setText(QString::number(Delta_E, 'f', 4));
    ui->param_Group1_p          ->setText(QString::number(p, 'f', 4));
    ui->param_Group1_Sigma_ap_d ->setText(QString::number(Sigma_ap_d, 'f', 4));
    ui->param_Group1_Sigma_ap_i ->setText(QString::number(Sigma_ap_i, 'f', 4));
}

void MainWindow::on__btn_CalcGroup_2_clicked()
{
    // Поршневая головка шатуна
    //  Расчет сечения А-А

    double Phi_msh, p_zdelta, p_0, F_p, Sigma_ap_a;
        Phi_msh     = ui->param_Group2_Phi_msh->text().toDouble();
        p_zdelta    = ui->param_Group2_p_zdelta->text().toDouble();
        p_0         = ui->param_Group2_p0->text().toDouble();
        F_p         = ui->param_Group2_Fp->text().toDouble();
        Sigma_ap_a  = ui->param_Group2_Sigma_ap_a->text().toDouble();

    double r_cp, d_n, m_p, lambda, Phi_shz, phi, k_sigma, eps_M, eps_N;

    double Omega, P_jp, R, N_j0, M_j0, N_j_phishz, M_j_phishz, K, F_g, F_v, Sigma_aj, P_sj, N_sj_phisg, M_sj_phisg,
            Sigma_acc, Sigma_max, Sigma_min, Sigma_m, Sigma_a, Sigma_ak;

    //  Расчет максимальной силы, растягивающей головку
    Omega = M_PI * n_N / 30;                       //n_N       ?
    r_cp = (d_g + d) / 4;
    F_g = (d_g + d) * l_sh;
    F_v = (d - d_n) * l_sh;
    K = E_sh * F_g / (E_sh * F_g + E_v * F_v);   //E_sh, E_v ?

        P_jp = -m_p * pow(Omega,2) * R * (1 + lambda);        //m_p, R, lambda    ?
        N_j0 = -P_jp * (0.572 - 0.0008 * Phi_shz);        //phi_shz           ?
        M_j0 = -P_jp * r_cp * (0.00033 * Phi_shz - 0.0297); //phi_shz           ?

        N_j_phishz = N_j0 * cos(Phi_shz) - 0.5 * P_jp * (sin(Phi_shz) - cos(Phi_shz));                       //phi_shz   ?
        M_j_phishz = M_j0 * r_cp * (1 - cos(Phi_shz) + 0.5 * P_jp * r_cp * (sin(Phi_shz) - cos(Phi_shz)));   //phi_shz   ?

        Sigma_aj = (2 * M_j_phishz * 6 * r_cp + h_g / (h_g * (2 * r_cp + h_g)) + K * N_j_phishz) * pow(10,-6) / (l_sh * h_g);

        P_sj = (p_zdelta - p_0) * F_p - m_p * R * pow(Omega,2) * (cos(phi) + lambda * cos(2 * phi)); //p_z_delta, p_0, F_p, m_p, R, phi, lambda ?

        //выбор слогаемых из таблиц
        N_sj_phisg = 3;
        M_sj_phisg = 7;
        // N_sj_phishz = P_sj * (N1 + (N2 - N3 - N4))
        // M_sj_phishz = P_sj * r_cp(M1 + (M2 - M3 - M4))

        Sigma_acc = (2 * M_sj_phisg * (6 * r_cp + h_g) / (h_g * (2 * r_cp + h_g)) + K * N_sj_phisg) * pow(10,-6) / (l_sh * h_g);
        Sigma_max = Sigma_ap_a + Sigma_aj;          //Sigma_ap_a ?
        Sigma_min = Sigma_ap_a + Sigma_acc;         //Sigma_ap_a ?
        Sigma_m = (Sigma_max + Sigma_min) / 2;
        Sigma_a = (Sigma_max - Sigma_min) / 2;
        Sigma_ak = Sigma_a * k_sigma / (eps_M * eps_N);

//    ui->param_Group2_
    ui->param_Group2_Omega      ->setText(QString::number(Omega, 'f', 4));
    ui->param_Group2_Pjp        ->setText(QString::number(P_jp, 'f', 4));
    ui->param_Group2_RSred      ->setText(QString::number(R, 'f', 4));
    ui->param_Group2_Nj0        ->setText(QString::number(N_j0, 'f', 4));
    ui->param_Group2_Mj0        ->setText(QString::number(M_j0, 'f', 4));
    ui->param_Group2_Nj         ->setText(QString::number(N_j_phishz, 'f', 4));
    ui->param_Group2_Mj         ->setText(QString::number(M_j_phishz, 'f', 4));
    ui->param_Group2_Nsj        ->setText(QString::number(N_sj_phisg, 'f', 4));
    ui->param_Group2_Msj        ->setText(QString::number(M_sj_phisg, 'f', 4));
    ui->param_Group2_K          ->setText(QString::number(K, 'f', 4));
    ui->param_Group2_Fg         ->setText(QString::number(F_g, 'f', 4));
    ui->param_Group2_Fv         ->setText(QString::number(F_v, 'f', 4));
    ui->param_Group2_Sigmaaj    ->setText(QString::number(Sigma_aj, 'f', 4));
    ui->param_Group2_Psg        ->setText(QString::number(P_sj, 'f', 4));
    ui->param_Group2_Sigma_acc  ->setText(QString::number(Sigma_acc, 'f', 4));
    ui->param_Group2_Sigma_max  ->setText(QString::number(Sigma_max, 'f', 4));
    ui->param_Group2_Sigma_min  ->setText(QString::number(Sigma_min, 'f', 4));
    ui->param_Group2_Sigma_m    ->setText(QString::number(Sigma_m, 'f', 4));
    ui->param_Group2_Sigma_a    ->setText(QString::number(Sigma_a, 'f', 4));
    ui->param_Group2_Sigma_ak   ->setText(QString::number(Sigma_ak, 'f', 4));
}

void MainWindow::on__btn_CalcGroup_3_clicked()
{
//    Кривошипная головка шатуна

    double m_sh, m_shp, m_shk;
    m_sh  = ui->param_Group3_msh->text().toDouble();
    m_shp = ui->param_Group3_mshp->text().toDouble();
    m_shk = ui->param_Group3_mshk->text().toDouble();

    double Omg_xx_max, R, m_p, lambda, Sigma_iz;

    double r_1, J_v, J, W_iz, F_g, m_kp, P_jp;

    r_1 = 0.5 * (d_2sh + 2 * t_v);
    J_v = l_k * pow(t_v,2);
    J = l_k * pow(0.5 * C_b - r_1,3);
    W_iz = l_k * pow(0.5 * C_b - r_1,2) / 6;
    F_g = 0.5 * l_k * (C_b - d_2sh);             // C_b   ?

    P_jp = pow(Omg_xx_max,2) * R * ((m_p + m_shp) * (1 + lambda) + (m_shk - m_kp)) * pow(10,-6);  // omg_xx_max, m_n, m_shp, lambda, m_shk, m_kp ?
    Sigma_iz = P_jp * (0.023 * C_b / ((1 + J_v / J) * W_iz) + 0.4 / F_g);

    ui->param_Group3_r1 ->setText(QString::number(r_1, 'f', 4));
    ui->param_Group3_Jv ->setText(QString::number(J_v, 'f', 4));
    ui->param_Group3_J  ->setText(QString::number(J, 'f', 4));
    ui->param_Group3_Wiz->setText(QString::number(W_iz, 'f', 4));
    ui->param_Group3_Fg ->setText(QString::number(F_g, 'f', 4));
    ui->param_Group3_mkr->setText(QString::number(m_kp, 'f', 4));
    ui->param_Group3_Pjp->setText(QString::number(P_jp, 'f', 4));
}

void MainWindow::on__btn_CalcDinam_clicked()
{
    double m_k, m_sh, m_p, R, omega, F_p;
    m_k     = ui->param_dinam_mk0   ->text().toDouble();
    m_sh    = ui->param_dinam_msht   ->text().toDouble();
    m_p     = ui->param_dinam_mp    ->text().toDouble();
    R       = ui->param_dinam_R     ->text().toDouble();
    omega   = ui->param_dinam_omega ->text().toDouble();
    F_p     = ui->param_dinam_Fp    ->text().toDouble();

    double m_shk, K_Rsh, K_Rk, K_Rsum;
    m_shk = 0.725 * m_sh;
    K_Rsh = -m_shk * R * pow(omega,2) / 1000;
    m_k   = m_k * F_p;
    K_Rk  = -m_k * R * pow(omega,2) / 1000;
    K_Rsum = K_Rsh + K_Rk;

    ui->param_dinam_mshk ->setText(QString::number(m_shk, 'f', 4));
    ui->param_dinam_KRsh ->setText(QString::number(K_Rsh, 'f', 4));
    ui->param_dinam_mk   ->setText(QString::number(m_k, 'f', 4));
    ui->param_dinam_KRk  ->setText(QString::number(K_Rk, 'f', 4));
    ui->param_dinam_KRsum->setText(QString::number(K_Rsum, 'f', 4));

    QStringList HeaderLabels;
    HeaderLabels << "φ(deg)"
                 << "Pi(φ)"
                 << "Pj(φ)"
                 << "PΣ(φ)";
    ui->table_dinam_Pphi->setColumnCount(4);
    ui->table_dinam_Pphi->setRowCount(73);
    ui->table_dinam_Pphi->setHorizontalHeaderLabels(HeaderLabels);

    QVector<double> Pi_var(74);
    int i=0;
    QFile file_Pi(":/Pi");
    if (!file_Pi.open(QIODevice::ReadOnly))
        return;

    QTextStream textStream_Pi(&file_Pi);

    while (!textStream_Pi.atEnd())
    {
        Pi_var[i] = textStream_Pi.readLine().toDouble();
        i++;
    }

    QVector<double> Pi(73), Pj(73), Psum(73);
    QVector<double> Deg(73);

    double m_post_, m_post, lambda = 0.1;
    m_post  = m_p + 1/3 * m_sh;
    m_post_ = m_post / F_p;

    double maxP = 0, minP = 0;

    for (int row=0; row<ui->table_dinam_Pphi->rowCount(); row++){
        double deg = M_PI / 180;
        double phi = row * 10 * deg;

        Pi[row] = 0;
        Pj[row] = 0;
        Psum[row] = 0;

        Deg[row] = row * 10;
//        double RND = rand() % 10000;
//        Pi[row] = RND;
        Pi[row] = Pi_var[row];
            if (Pi[row] > maxP) maxP = Pi[row];
            if (Pi[row] < minP) minP = Pi[row];
        Pj[row] = m_post_ * R * pow(omega,1) * (cos(Deg[row]) + lambda * cos(2 * Deg[row]));
        Pj[row] = Pj[row]/10000;
            if (Pj[row] > maxP) maxP = Pj[row];
            if (Pj[row] < minP) minP = Pj[row];
        Psum[row] = Pj[row] + Pi[row];
            if (Psum[row] > maxP) maxP = Psum[row];
            if (Psum[row] < minP) minP = Psum[row];

        for (int col=0; col<ui->table_dinam_Pphi->columnCount(); col++){
            QTableWidgetItem *Item = new QTableWidgetItem();
            double VarItem;
            switch (col) {
            case 0:
                VarItem = Deg[row];
                break;
            case 1:
                VarItem = Pi[row];
                break;
            case 2:
                VarItem = Pj[row];
                break;
            case 3:
                VarItem = Psum[row];
                break;
            default:
                break;
            }
            if (col == 0)
                Item->setText(QString::number((int)VarItem));
            else
                Item->setText(QString::number(VarItem, 'f', 4));
                Item->setTextAlignment(Qt::AlignCenter);
            ui->table_dinam_Pphi->setItem(row,col,Item);
            ui->table_dinam_Pphi->setColumnWidth(col,78);
        }
    }

    ui->plot_dinam_Pphi->clearGraphs();
    ui->plot_dinam_Pphi->addGraph();
    ui->plot_dinam_Pphi->graph(0)->setData(Deg,Pj);
    ui->plot_dinam_Pphi->graph(0)->setLineStyle((QCPGraph::lsLine));
    ui->plot_dinam_Pphi->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot_dinam_Pphi->graph(0)->setPen(QPen(Qt::black));
    ui->plot_dinam_Pphi->graph(0)->setName("Pj");

    ui->plot_dinam_Pphi->addGraph();
    ui->plot_dinam_Pphi->graph(1)->setData(Deg,Pi);
    ui->plot_dinam_Pphi->graph(1)->setLineStyle((QCPGraph::lsLine));
    ui->plot_dinam_Pphi->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot_dinam_Pphi->graph(1)->setPen(QPen(Qt::blue));
    ui->plot_dinam_Pphi->graph(1)->setName("Pi");

    ui->plot_dinam_Pphi->addGraph();
    ui->plot_dinam_Pphi->graph(2)->setData(Deg,Psum);
    ui->plot_dinam_Pphi->graph(2)->setLineStyle((QCPGraph::lsLine));
    ui->plot_dinam_Pphi->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot_dinam_Pphi->graph(2)->setPen(QPen(Qt::red));
    ui->plot_dinam_Pphi->graph(2)->setName("Psum");

    ui->plot_dinam_Pphi->xAxis->setLabel("φ");
    ui->plot_dinam_Pphi->xAxis->setRange(0, 720);
    ui->plot_dinam_Pphi->yAxis->setRange(1.05*minP, 1.05*maxP);
    ui->plot_dinam_Pphi->replot();

    ui->plot_dinam_Pphi->legend->setVisible(true);
    ui->plot_dinam_Pphi->legend->setBrush(QBrush(QColor(255,255,255,150)));
    ui->plot_dinam_Pphi->yAxis->grid()->setSubGridVisible(true);
    ui->plot_dinam_Pphi->xAxis->grid()->setSubGridVisible(true);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    QVector<double> d1(74), d2(74), d3(74), d4(74);
    i=0;
    QFile file_d1(":/d1"),
            file_d2(":/d2"),
            file_d3(":/d3"),
            file_d4(":/d4");
    if (!file_d1.open(QIODevice::ReadOnly) ||
            !file_d2.open(QIODevice::ReadOnly) ||
            !file_d3.open(QIODevice::ReadOnly) ||
            !file_d4.open(QIODevice::ReadOnly))
        return;

    QTextStream textStream_d1(&file_d1);
    QTextStream textStream_d2(&file_d2);
    QTextStream textStream_d3(&file_d3);
    QTextStream textStream_d4(&file_d4);

    while (!textStream_d1.atEnd() ||
           !textStream_d2.atEnd() ||
           !textStream_d3.atEnd() ||
           !textStream_d4.atEnd())
    {
        d1[i] = textStream_d1.readLine().toDouble();
        d2[i] = textStream_d2.readLine().toDouble();
        d3[i] = textStream_d3.readLine().toDouble();
        d4[i] = textStream_d4.readLine().toDouble();
        i++;
    }

    QVector<double> N(73), S(73), K(73), T(73);
    double Nmin=0, Nmax=0;
    double Smin=0, Smax=0;
    double Kmin=0, Kmax=0;
    double Tmin=0, Tmax=0;

    for (int i = 0; i < 73; i++) {
        N[i] = Psum[i] * d1[i];
        S[i] = Psum[i] * d2[i];
        K[i] = Psum[i] * d3[i];
        T[i] = Psum[i] * d4[i];

        if (N[i] > Nmax) Nmax = N[i];
        if (N[i] < Nmin) Nmin = N[i];

        if (S[i] > Smax) Smax = S[i];
        if (S[i] < Smin) Smin = S[i];

        if (K[i] > Kmax) Kmax = K[i];
        if (K[i] < Kmin) Kmin = K[i];

        if (T[i] > Tmax) Tmax = T[i];
        if (T[i] < Tmin) Tmin = T[i];
    }

    ui->plot_dinam_N->clearGraphs();
    ui->plot_dinam_N->addGraph();
    ui->plot_dinam_N->graph(0)->setData(Deg,N);
    ui->plot_dinam_N->graph(0)->setLineStyle((QCPGraph::lsLine));
//    ui->plot_dinam_N->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot_dinam_N->xAxis->setLabel("φ");
    ui->plot_dinam_N->yAxis->setLabel("N(φ)");
    ui->plot_dinam_N->xAxis->setRange(0, 720);
    ui->plot_dinam_N->yAxis->setRange(Nmin - 0.05 * Nmin, Nmax + 0.05 * Nmax);
    ui->plot_dinam_N->replot();

    ui->plot_dinam_S->clearGraphs();
    ui->plot_dinam_S->addGraph();
    ui->plot_dinam_S->graph(0)->setData(Deg,S);
    ui->plot_dinam_S->graph(0)->setLineStyle((QCPGraph::lsLine));
//    ui->plot_dinam_S->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot_dinam_S->xAxis->setLabel("φ");
    ui->plot_dinam_S->yAxis->setLabel("S(φ)");
    ui->plot_dinam_S->xAxis->setRange(0, 720);
    ui->plot_dinam_S->yAxis->setRange(Smin - 0.05 * Smin, Smax + 0.05 * Smax);
    ui->plot_dinam_S->replot();

    ui->plot_dinam_K->clearGraphs();
    ui->plot_dinam_K->addGraph();
    ui->plot_dinam_K->graph(0)->setData(Deg,K);
    ui->plot_dinam_K->graph(0)->setLineStyle((QCPGraph::lsLine));
//    ui->plot_dinam_K->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot_dinam_K->xAxis->setLabel("φ");
    ui->plot_dinam_K->yAxis->setLabel("K(φ)");
    ui->plot_dinam_K->xAxis->setRange(0, 720);
    ui->plot_dinam_K->yAxis->setRange(Kmin - 0.05 * Kmin, Kmax + 0.05 * Kmax);
    ui->plot_dinam_K->replot();

    ui->plot_dinam_T->clearGraphs();
    ui->plot_dinam_T->addGraph();
    ui->plot_dinam_T->graph(0)->setData(Deg,T);
    ui->plot_dinam_T->graph(0)->setLineStyle(QCPGraph::lsLine);
//    ui->plot_dinam_T->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->plot_dinam_T->xAxis->setLabel("φ");
    ui->plot_dinam_T->yAxis->setLabel("T(φ)");
    ui->plot_dinam_T->xAxis->setRange(0, 720);
    ui->plot_dinam_T->yAxis->setRange(Tmin - 0.05 * Tmin, Tmax + 0.05 * Tmax);
    ui->plot_dinam_T->replot();


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    
}
