OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8795348) q[0];
sx q[0];
rz(-1.4095925) q[0];
sx q[0];
rz(-1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(6.8847818) q[1];
sx q[1];
rz(9.8431982) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19791662) q[0];
sx q[0];
rz(-2.0352053) q[0];
sx q[0];
rz(2.8465413) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7069874) q[2];
sx q[2];
rz(-1.681466) q[2];
sx q[2];
rz(1.9905123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99416713) q[1];
sx q[1];
rz(-0.52417437) q[1];
sx q[1];
rz(-2.0736573) q[1];
x q[2];
rz(-1.1316142) q[3];
sx q[3];
rz(-1.3703128) q[3];
sx q[3];
rz(2.2905614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-2.3036172) q[2];
rz(0.49301246) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-0.078991927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.2778506) q[0];
rz(0.17678075) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(0.4321672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55789253) q[0];
sx q[0];
rz(-1.0520792) q[0];
sx q[0];
rz(1.2516663) q[0];
rz(-0.400153) q[2];
sx q[2];
rz(-1.1740985) q[2];
sx q[2];
rz(-2.951159) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13465263) q[1];
sx q[1];
rz(-1.4538527) q[1];
sx q[1];
rz(-2.0140531) q[1];
rz(-1.5315941) q[3];
sx q[3];
rz(-0.56534846) q[3];
sx q[3];
rz(2.431734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5923578) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(0.48669997) q[2];
rz(1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(-2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040722672) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(1.4886645) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(2.6352077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8538118) q[0];
sx q[0];
rz(-1.7475442) q[0];
sx q[0];
rz(1.9275097) q[0];
rz(-pi) q[1];
rz(1.330553) q[2];
sx q[2];
rz(-2.2908387) q[2];
sx q[2];
rz(2.4400997) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80458927) q[1];
sx q[1];
rz(-0.5491921) q[1];
sx q[1];
rz(-0.6188436) q[1];
rz(-pi) q[2];
rz(0.18249986) q[3];
sx q[3];
rz(-1.9347408) q[3];
sx q[3];
rz(-2.5516627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69904077) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(0.59147269) q[2];
rz(-2.5555723) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(-1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(0.83918321) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(1.5485839) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73491053) q[0];
sx q[0];
rz(-1.6618177) q[0];
sx q[0];
rz(-1.1199561) q[0];
rz(-pi) q[1];
rz(2.2976539) q[2];
sx q[2];
rz(-0.91274777) q[2];
sx q[2];
rz(2.7968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81739391) q[1];
sx q[1];
rz(-0.84016582) q[1];
sx q[1];
rz(-2.6170931) q[1];
rz(-0.56460103) q[3];
sx q[3];
rz(-2.6435916) q[3];
sx q[3];
rz(2.6548487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-3.0419066) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(-1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(2.8856522) q[0];
rz(0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(-0.76006132) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81239031) q[0];
sx q[0];
rz(-0.15604067) q[0];
sx q[0];
rz(2.4183256) q[0];
x q[1];
rz(-3.0779482) q[2];
sx q[2];
rz(-1.1568937) q[2];
sx q[2];
rz(-2.6457583) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0415503) q[1];
sx q[1];
rz(-1.6987545) q[1];
sx q[1];
rz(0.12771878) q[1];
rz(-0.92442583) q[3];
sx q[3];
rz(-1.8674208) q[3];
sx q[3];
rz(-2.7588206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57050675) q[2];
sx q[2];
rz(-1.3723624) q[2];
sx q[2];
rz(-0.6742397) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6102585) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(-1.1791139) q[0];
rz(0.20482652) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(-2.0746322) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5087591) q[0];
sx q[0];
rz(-1.5501223) q[0];
sx q[0];
rz(1.585292) q[0];
x q[1];
rz(2.009379) q[2];
sx q[2];
rz(-1.1819934) q[2];
sx q[2];
rz(0.68225451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9251717) q[1];
sx q[1];
rz(-0.74837084) q[1];
sx q[1];
rz(-1.5009576) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1596998) q[3];
sx q[3];
rz(-0.80280639) q[3];
sx q[3];
rz(-1.6646977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.431488) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(2.4152749) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(-0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2840435) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(-1.7204826) q[0];
rz(2.9395318) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(-0.85817671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1976397) q[0];
sx q[0];
rz(-1.4240992) q[0];
sx q[0];
rz(-3.0118045) q[0];
x q[1];
rz(1.6740587) q[2];
sx q[2];
rz(-2.729136) q[2];
sx q[2];
rz(2.7445284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9550025) q[1];
sx q[1];
rz(-1.5860671) q[1];
sx q[1];
rz(-3.112622) q[1];
x q[2];
rz(1.2248366) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(-2.791415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28785607) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637852) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(-2.7668787) q[0];
rz(2.162714) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(-1.3495061) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9106306) q[0];
sx q[0];
rz(-2.1719143) q[0];
sx q[0];
rz(-2.0448951) q[0];
rz(-1.5937514) q[2];
sx q[2];
rz(-2.288753) q[2];
sx q[2];
rz(-0.1644451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.046557758) q[1];
sx q[1];
rz(-1.7525502) q[1];
sx q[1];
rz(3.0912116) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3091062) q[3];
sx q[3];
rz(-1.040254) q[3];
sx q[3];
rz(-1.7390651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4902041) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(-1.6040241) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8309098) q[0];
sx q[0];
rz(-1.3498107) q[0];
sx q[0];
rz(1.2587147) q[0];
rz(-pi) q[1];
rz(-1.6749009) q[2];
sx q[2];
rz(-1.1541919) q[2];
sx q[2];
rz(2.5660851) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4496085) q[1];
sx q[1];
rz(-2.2712628) q[1];
sx q[1];
rz(-1.4222172) q[1];
x q[2];
rz(-2.2213307) q[3];
sx q[3];
rz(-0.55763054) q[3];
sx q[3];
rz(-2.695431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1016772) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-0.76134479) q[2];
rz(-2.2411761) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3392357) q[0];
sx q[0];
rz(-0.46827066) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(0.63198173) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8666329) q[0];
sx q[0];
rz(-1.8909591) q[0];
sx q[0];
rz(0.76036705) q[0];
rz(-pi) q[1];
rz(-1.1502613) q[2];
sx q[2];
rz(-1.0610126) q[2];
sx q[2];
rz(-1.2812986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6001544) q[1];
sx q[1];
rz(-0.46240515) q[1];
sx q[1];
rz(0.17141436) q[1];
x q[2];
rz(-0.93654376) q[3];
sx q[3];
rz(-1.9359971) q[3];
sx q[3];
rz(1.4432424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0163991) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(0.94474244) q[2];
rz(-0.38481209) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(-0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(0.8846994) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-2.5169218) q[2];
sx q[2];
rz(-2.2041337) q[2];
sx q[2];
rz(0.49908257) q[2];
rz(1.9839722) q[3];
sx q[3];
rz(-0.4784085) q[3];
sx q[3];
rz(2.9984409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
