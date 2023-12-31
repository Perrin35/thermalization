OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6611377) q[0];
sx q[0];
rz(-1.2167131) q[0];
sx q[0];
rz(0.35868355) q[0];
rz(-1.0010927) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(0.18730883) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5536082) q[1];
sx q[1];
rz(-1.6349287) q[1];
sx q[1];
rz(1.1448121) q[1];
x q[2];
rz(2.7798153) q[3];
sx q[3];
rz(-0.28218111) q[3];
sx q[3];
rz(0.85229814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-3.1100173) q[2];
rz(-1.8850373) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(-2.6020715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8115494) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(-0.051785843) q[0];
x q[1];
rz(1.9794481) q[2];
sx q[2];
rz(-0.89102972) q[2];
sx q[2];
rz(-1.0811999) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.19903781) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(-1.1953137) q[1];
rz(-0.86124729) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(1.4351821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3617525) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(-1.386863) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(-0.2562491) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1218131) q[0];
sx q[0];
rz(-2.461713) q[0];
sx q[0];
rz(2.7133184) q[0];
rz(0.44552866) q[2];
sx q[2];
rz(-1.1330714) q[2];
sx q[2];
rz(-2.90403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5970522) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(2.271133) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6085298) q[3];
sx q[3];
rz(-1.3940666) q[3];
sx q[3];
rz(-2.7023774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(1.3556708) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(0.25948778) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(0.73192275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67128348) q[0];
sx q[0];
rz(-1.0571612) q[0];
sx q[0];
rz(2.1156103) q[0];
x q[1];
rz(-0.71728431) q[2];
sx q[2];
rz(-1.4826164) q[2];
sx q[2];
rz(-0.35536534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6432453) q[1];
sx q[1];
rz(-2.8162662) q[1];
sx q[1];
rz(-2.494032) q[1];
rz(-0.51283522) q[3];
sx q[3];
rz(-2.8898015) q[3];
sx q[3];
rz(-2.4076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(1.2623825) q[0];
rz(1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-2.343822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.939643) q[0];
sx q[0];
rz(-1.5683096) q[0];
sx q[0];
rz(1.6909088) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8380978) q[2];
sx q[2];
rz(-1.7804838) q[2];
sx q[2];
rz(1.4022624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7792015) q[1];
sx q[1];
rz(-2.5070842) q[1];
sx q[1];
rz(-1.7290551) q[1];
rz(0.4425211) q[3];
sx q[3];
rz(-1.2708775) q[3];
sx q[3];
rz(2.4181441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(-0.30203715) q[2];
rz(-1.9942412) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.9858032) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-0.070080431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2698343) q[0];
sx q[0];
rz(-1.4176798) q[0];
sx q[0];
rz(-1.7437177) q[0];
rz(2.5885133) q[2];
sx q[2];
rz(-2.2773909) q[2];
sx q[2];
rz(1.4097708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8225704) q[1];
sx q[1];
rz(-2.550107) q[1];
sx q[1];
rz(0.2925847) q[1];
rz(3.1061884) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30248102) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(-0.85987464) q[0];
rz(-1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(-0.0079356114) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82848362) q[0];
sx q[0];
rz(-2.5223753) q[0];
sx q[0];
rz(-2.6909268) q[0];
rz(-pi) q[1];
rz(-2.7107312) q[2];
sx q[2];
rz(-2.2164946) q[2];
sx q[2];
rz(0.25260392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7759526) q[1];
sx q[1];
rz(-1.0675634) q[1];
sx q[1];
rz(2.3535437) q[1];
rz(-pi) q[2];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.946297) q[3];
sx q[3];
rz(0.34208959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(-2.725214) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96173441) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(-0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.6392802) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4541894) q[0];
sx q[0];
rz(-0.25035509) q[0];
sx q[0];
rz(3.1012015) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6697568) q[2];
sx q[2];
rz(-2.2484357) q[2];
sx q[2];
rz(-0.0080646947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8870526) q[1];
sx q[1];
rz(-2.2215543) q[1];
sx q[1];
rz(-0.73144967) q[1];
rz(-pi) q[2];
rz(-1.1542529) q[3];
sx q[3];
rz(-0.22437469) q[3];
sx q[3];
rz(2.8191872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97312462) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(2.705943) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(2.7246144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57396736) q[0];
sx q[0];
rz(-1.5331393) q[0];
sx q[0];
rz(-2.2249939) q[0];
rz(-pi) q[1];
rz(-0.68848227) q[2];
sx q[2];
rz(-2.3496029) q[2];
sx q[2];
rz(-0.85306963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3560564) q[1];
sx q[1];
rz(-1.0724663) q[1];
sx q[1];
rz(2.9836125) q[1];
x q[2];
rz(-2.3973893) q[3];
sx q[3];
rz(-2.9508698) q[3];
sx q[3];
rz(2.423559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(0.22582516) q[2];
rz(0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(1.5243994) q[0];
rz(-0.95364755) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1117301) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(-1.8882621) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1216713) q[2];
sx q[2];
rz(-1.7977062) q[2];
sx q[2];
rz(-2.801148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9433371) q[1];
sx q[1];
rz(-0.88611929) q[1];
sx q[1];
rz(-0.68590045) q[1];
rz(-pi) q[2];
rz(0.15354746) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(-0.49680199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9983457) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(1.5608816) q[2];
sx q[2];
rz(-2.0799939) q[2];
sx q[2];
rz(-1.6891198) q[2];
rz(0.090311269) q[3];
sx q[3];
rz(-1.3406546) q[3];
sx q[3];
rz(-1.8484074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
