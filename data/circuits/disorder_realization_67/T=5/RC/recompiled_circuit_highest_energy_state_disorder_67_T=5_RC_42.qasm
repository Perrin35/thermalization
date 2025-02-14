OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8494868) q[0];
sx q[0];
rz(-2.5295244) q[0];
sx q[0];
rz(3.0906313) q[0];
rz(-1.851097) q[1];
sx q[1];
rz(4.4176997) q[1];
sx q[1];
rz(7.1292276) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.030662) q[0];
sx q[0];
rz(-1.7314859) q[0];
sx q[0];
rz(2.7881505) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92535352) q[2];
sx q[2];
rz(-0.45956372) q[2];
sx q[2];
rz(-1.1499576) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4197246) q[1];
sx q[1];
rz(-1.8375735) q[1];
sx q[1];
rz(0.50290033) q[1];
rz(-pi) q[2];
rz(1.5624763) q[3];
sx q[3];
rz(-1.6826311) q[3];
sx q[3];
rz(-0.40598265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35752615) q[2];
sx q[2];
rz(-2.1937723) q[2];
sx q[2];
rz(2.2785462) q[2];
rz(1.1132025) q[3];
sx q[3];
rz(-2.2200255) q[3];
sx q[3];
rz(-0.79871261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5384101) q[0];
sx q[0];
rz(-0.69284678) q[0];
sx q[0];
rz(-2.8811654) q[0];
rz(0.32568359) q[1];
sx q[1];
rz(-2.1574056) q[1];
sx q[1];
rz(-0.30776417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9887675) q[0];
sx q[0];
rz(-0.77251311) q[0];
sx q[0];
rz(1.9475742) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57450104) q[2];
sx q[2];
rz(-2.5972453) q[2];
sx q[2];
rz(-0.33657227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.766749) q[1];
sx q[1];
rz(-1.6628712) q[1];
sx q[1];
rz(-1.8122775) q[1];
rz(-0.006853718) q[3];
sx q[3];
rz(-2.2061976) q[3];
sx q[3];
rz(-1.3069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3088358) q[2];
sx q[2];
rz(-0.86897659) q[2];
sx q[2];
rz(-2.5858509) q[2];
rz(-0.51221383) q[3];
sx q[3];
rz(-0.54849505) q[3];
sx q[3];
rz(0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.888805) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(-2.3922701) q[0];
rz(-1.0961756) q[1];
sx q[1];
rz(-1.3744033) q[1];
sx q[1];
rz(-0.38415092) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6644728) q[0];
sx q[0];
rz(-0.64711231) q[0];
sx q[0];
rz(0.44575925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60127705) q[2];
sx q[2];
rz(-1.3601175) q[2];
sx q[2];
rz(0.36128584) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.322142) q[1];
sx q[1];
rz(-1.8729728) q[1];
sx q[1];
rz(-3.0858529) q[1];
x q[2];
rz(-0.29112737) q[3];
sx q[3];
rz(-2.2343221) q[3];
sx q[3];
rz(-0.71602589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38203794) q[2];
sx q[2];
rz(-1.8363154) q[2];
sx q[2];
rz(-1.6273512) q[2];
rz(-2.4115327) q[3];
sx q[3];
rz(-2.4000013) q[3];
sx q[3];
rz(-0.39456427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6877614) q[0];
sx q[0];
rz(-1.4241968) q[0];
sx q[0];
rz(0.060039595) q[0];
rz(-0.91375786) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(-2.2027016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.765678) q[0];
sx q[0];
rz(-1.5374827) q[0];
sx q[0];
rz(0.062996431) q[0];
rz(-pi) q[1];
rz(1.6020158) q[2];
sx q[2];
rz(-1.6531715) q[2];
sx q[2];
rz(2.9766469) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9524433) q[1];
sx q[1];
rz(-1.5464142) q[1];
sx q[1];
rz(-1.0827176) q[1];
rz(-pi) q[2];
rz(-2.4603953) q[3];
sx q[3];
rz(-1.0170611) q[3];
sx q[3];
rz(-0.21986976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1303723) q[2];
sx q[2];
rz(-0.29272407) q[2];
sx q[2];
rz(2.4893153) q[2];
rz(-1.8782714) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(-1.9704845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0703099) q[0];
sx q[0];
rz(-1.2391397) q[0];
sx q[0];
rz(-0.5624482) q[0];
rz(-1.3485472) q[1];
sx q[1];
rz(-2.8548073) q[1];
sx q[1];
rz(-2.5602692) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59636694) q[0];
sx q[0];
rz(-2.2196182) q[0];
sx q[0];
rz(-2.5741386) q[0];
x q[1];
rz(1.6069218) q[2];
sx q[2];
rz(-1.110552) q[2];
sx q[2];
rz(2.4789916) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5229458) q[1];
sx q[1];
rz(-0.60877555) q[1];
sx q[1];
rz(2.2227915) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64251803) q[3];
sx q[3];
rz(-1.5356365) q[3];
sx q[3];
rz(-0.25504756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44744667) q[2];
sx q[2];
rz(-1.7116825) q[2];
sx q[2];
rz(-0.68667975) q[2];
rz(-0.40003362) q[3];
sx q[3];
rz(-1.2587849) q[3];
sx q[3];
rz(0.011628477) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069063) q[0];
sx q[0];
rz(-2.1676368) q[0];
sx q[0];
rz(0.45027012) q[0];
rz(0.88092583) q[1];
sx q[1];
rz(-1.9073146) q[1];
sx q[1];
rz(-1.1164104) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085297708) q[0];
sx q[0];
rz(-1.3816815) q[0];
sx q[0];
rz(-2.3448639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47412761) q[2];
sx q[2];
rz(-2.7517786) q[2];
sx q[2];
rz(1.4963983) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7959545) q[1];
sx q[1];
rz(-1.5511723) q[1];
sx q[1];
rz(2.1810637) q[1];
rz(-pi) q[2];
rz(0.34120543) q[3];
sx q[3];
rz(-1.998105) q[3];
sx q[3];
rz(-0.36580929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66699666) q[2];
sx q[2];
rz(-0.23385364) q[2];
sx q[2];
rz(-2.4901566) q[2];
rz(2.3598119) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(-0.93530161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5485789) q[0];
sx q[0];
rz(-2.9055556) q[0];
sx q[0];
rz(2.6759942) q[0];
rz(-1.9904526) q[1];
sx q[1];
rz(-1.2460037) q[1];
sx q[1];
rz(1.8021072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5975153) q[0];
sx q[0];
rz(-1.9751514) q[0];
sx q[0];
rz(-2.5687575) q[0];
x q[1];
rz(-1.5268099) q[2];
sx q[2];
rz(-1.5325452) q[2];
sx q[2];
rz(1.418131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7581343) q[1];
sx q[1];
rz(-1.3515359) q[1];
sx q[1];
rz(3.0644528) q[1];
rz(0.42134704) q[3];
sx q[3];
rz(-0.94910062) q[3];
sx q[3];
rz(-2.5558215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0307978) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(0.01290713) q[2];
rz(-3.0068093) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(0.24136647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6808692) q[0];
sx q[0];
rz(-0.091982059) q[0];
sx q[0];
rz(-1.9665834) q[0];
rz(2.3911632) q[1];
sx q[1];
rz(-1.78396) q[1];
sx q[1];
rz(-2.7480385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1712883) q[0];
sx q[0];
rz(-1.7455202) q[0];
sx q[0];
rz(-2.1615684) q[0];
rz(-pi) q[1];
rz(2.4855087) q[2];
sx q[2];
rz(-1.4781845) q[2];
sx q[2];
rz(-0.54457918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8361349) q[1];
sx q[1];
rz(-1.4609481) q[1];
sx q[1];
rz(-2.0019462) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2116488) q[3];
sx q[3];
rz(-2.5287147) q[3];
sx q[3];
rz(1.5811435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5081818) q[2];
sx q[2];
rz(-1.1339302) q[2];
sx q[2];
rz(0.22937648) q[2];
rz(0.091382787) q[3];
sx q[3];
rz(-1.4437851) q[3];
sx q[3];
rz(-2.5858333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2185739) q[0];
sx q[0];
rz(-0.78892437) q[0];
sx q[0];
rz(-2.5564585) q[0];
rz(1.2773889) q[1];
sx q[1];
rz(-1.9265415) q[1];
sx q[1];
rz(0.1543943) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.383787) q[0];
sx q[0];
rz(-0.4497779) q[0];
sx q[0];
rz(-2.8881254) q[0];
rz(-pi) q[1];
rz(-1.2227258) q[2];
sx q[2];
rz(-1.1670975) q[2];
sx q[2];
rz(-2.8793546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33340633) q[1];
sx q[1];
rz(-1.8318614) q[1];
sx q[1];
rz(-0.33269791) q[1];
rz(-pi) q[2];
rz(-2.0386857) q[3];
sx q[3];
rz(-1.4135191) q[3];
sx q[3];
rz(-3.0426971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1207017) q[2];
sx q[2];
rz(-0.88940826) q[2];
sx q[2];
rz(-2.6013539) q[2];
rz(-2.7464416) q[3];
sx q[3];
rz(-0.65427798) q[3];
sx q[3];
rz(-1.3490217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9304792) q[0];
sx q[0];
rz(-1.2808639) q[0];
sx q[0];
rz(-1.6092009) q[0];
rz(-0.92597517) q[1];
sx q[1];
rz(-1.7151058) q[1];
sx q[1];
rz(-1.4769953) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4173399) q[0];
sx q[0];
rz(-0.84463929) q[0];
sx q[0];
rz(-1.5717616) q[0];
rz(-pi) q[1];
rz(2.8468644) q[2];
sx q[2];
rz(-1.9325614) q[2];
sx q[2];
rz(-0.28361646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.033869318) q[1];
sx q[1];
rz(-2.1185912) q[1];
sx q[1];
rz(-2.7366927) q[1];
rz(-0.98201237) q[3];
sx q[3];
rz(-1.9470306) q[3];
sx q[3];
rz(0.8134977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28107873) q[2];
sx q[2];
rz(-1.8033359) q[2];
sx q[2];
rz(2.4930084) q[2];
rz(2.4159238) q[3];
sx q[3];
rz(-0.23894335) q[3];
sx q[3];
rz(1.749595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3758748) q[0];
sx q[0];
rz(-0.97828843) q[0];
sx q[0];
rz(0.82213415) q[0];
rz(-1.0646461) q[1];
sx q[1];
rz(-1.4551661) q[1];
sx q[1];
rz(-1.0600769) q[1];
rz(-2.9072472) q[2];
sx q[2];
rz(-0.6498944) q[2];
sx q[2];
rz(1.9882974) q[2];
rz(1.0242994) q[3];
sx q[3];
rz(-2.306675) q[3];
sx q[3];
rz(-2.6425895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
