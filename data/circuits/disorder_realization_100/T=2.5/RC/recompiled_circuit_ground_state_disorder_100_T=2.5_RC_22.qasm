OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.36432) q[0];
sx q[0];
rz(3.9914208) q[0];
sx q[0];
rz(10.963647) q[0];
rz(-0.76378769) q[1];
sx q[1];
rz(3.0883767) q[1];
sx q[1];
rz(8.5903066) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0758103) q[0];
sx q[0];
rz(-2.6666323) q[0];
sx q[0];
rz(1.6823099) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1727512) q[2];
sx q[2];
rz(-1.0263832) q[2];
sx q[2];
rz(2.9679342) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.149772) q[1];
sx q[1];
rz(-1.705994) q[1];
sx q[1];
rz(2.1674431) q[1];
rz(-pi) q[2];
rz(2.8230366) q[3];
sx q[3];
rz(-1.1906822) q[3];
sx q[3];
rz(1.8442228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46204627) q[2];
sx q[2];
rz(-1.5507689) q[2];
sx q[2];
rz(-2.4756883) q[2];
rz(2.7146961) q[3];
sx q[3];
rz(-0.96564966) q[3];
sx q[3];
rz(0.20195937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17215984) q[0];
sx q[0];
rz(-0.96123022) q[0];
sx q[0];
rz(-0.52312553) q[0];
rz(-0.41921774) q[1];
sx q[1];
rz(-2.9328177) q[1];
sx q[1];
rz(2.8255393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8652865) q[0];
sx q[0];
rz(-1.3973736) q[0];
sx q[0];
rz(0.68626173) q[0];
x q[1];
rz(-2.9954073) q[2];
sx q[2];
rz(-1.5954895) q[2];
sx q[2];
rz(0.05505148) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16686116) q[1];
sx q[1];
rz(-0.63878794) q[1];
sx q[1];
rz(2.0898934) q[1];
rz(-pi) q[2];
rz(2.7160268) q[3];
sx q[3];
rz(-1.4168973) q[3];
sx q[3];
rz(2.9776412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54059625) q[2];
sx q[2];
rz(-2.5935846) q[2];
sx q[2];
rz(2.0151286) q[2];
rz(-0.90538853) q[3];
sx q[3];
rz(-2.0784046) q[3];
sx q[3];
rz(-1.6006914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044540502) q[0];
sx q[0];
rz(-1.4690761) q[0];
sx q[0];
rz(-1.4918406) q[0];
rz(2.7482694) q[1];
sx q[1];
rz(-1.2231188) q[1];
sx q[1];
rz(0.16600674) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9764023) q[0];
sx q[0];
rz(-0.87346625) q[0];
sx q[0];
rz(1.1572474) q[0];
rz(-pi) q[1];
rz(2.2691378) q[2];
sx q[2];
rz(-1.7042245) q[2];
sx q[2];
rz(-0.39233559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17540784) q[1];
sx q[1];
rz(-1.3864751) q[1];
sx q[1];
rz(2.4173474) q[1];
rz(-pi) q[2];
rz(0.29434926) q[3];
sx q[3];
rz(-2.6913319) q[3];
sx q[3];
rz(-0.0014622268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7321221) q[2];
sx q[2];
rz(-0.25264838) q[2];
sx q[2];
rz(0.18969336) q[2];
rz(2.7267552) q[3];
sx q[3];
rz(-1.8388351) q[3];
sx q[3];
rz(0.10492575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6787427) q[0];
sx q[0];
rz(-2.3609556) q[0];
sx q[0];
rz(-1.8829874) q[0];
rz(-1.3858093) q[1];
sx q[1];
rz(-2.7744727) q[1];
sx q[1];
rz(2.3585336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.258062) q[0];
sx q[0];
rz(-1.6029164) q[0];
sx q[0];
rz(1.5954992) q[0];
rz(-pi) q[1];
x q[1];
rz(0.003764228) q[2];
sx q[2];
rz(-2.3329511) q[2];
sx q[2];
rz(0.54169858) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2338581) q[1];
sx q[1];
rz(-2.0480262) q[1];
sx q[1];
rz(1.0027755) q[1];
rz(-3.0659292) q[3];
sx q[3];
rz(-0.35705413) q[3];
sx q[3];
rz(0.018659534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9328402) q[2];
sx q[2];
rz(-1.9806661) q[2];
sx q[2];
rz(-1.0932385) q[2];
rz(0.42539635) q[3];
sx q[3];
rz(-1.2109816) q[3];
sx q[3];
rz(2.0591327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0699186) q[0];
sx q[0];
rz(-2.2921083) q[0];
sx q[0];
rz(1.3776394) q[0];
rz(-2.2397974) q[1];
sx q[1];
rz(-1.8865562) q[1];
sx q[1];
rz(2.3756383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9760343) q[0];
sx q[0];
rz(-1.9355825) q[0];
sx q[0];
rz(2.4976612) q[0];
rz(-2.7904205) q[2];
sx q[2];
rz(-0.72423705) q[2];
sx q[2];
rz(1.8623095) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88845612) q[1];
sx q[1];
rz(-0.52475861) q[1];
sx q[1];
rz(-2.8349459) q[1];
x q[2];
rz(-0.28217989) q[3];
sx q[3];
rz(-1.9934142) q[3];
sx q[3];
rz(-2.9423892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19096601) q[2];
sx q[2];
rz(-2.2910255) q[2];
sx q[2];
rz(-1.4858049) q[2];
rz(-1.1367669) q[3];
sx q[3];
rz(-0.4509238) q[3];
sx q[3];
rz(0.76751149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9039827) q[0];
sx q[0];
rz(-0.5466277) q[0];
sx q[0];
rz(0.23608635) q[0];
rz(-1.4769332) q[1];
sx q[1];
rz(-1.6496941) q[1];
sx q[1];
rz(1.9794855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8315036) q[0];
sx q[0];
rz(-0.79206249) q[0];
sx q[0];
rz(-2.8312324) q[0];
rz(1.4496153) q[2];
sx q[2];
rz(-1.9773537) q[2];
sx q[2];
rz(0.32025619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9368105) q[1];
sx q[1];
rz(-0.46012649) q[1];
sx q[1];
rz(-2.2304753) q[1];
rz(1.5959019) q[3];
sx q[3];
rz(-1.5179765) q[3];
sx q[3];
rz(-2.3642295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7621925) q[2];
sx q[2];
rz(-0.47319943) q[2];
sx q[2];
rz(3.1326478) q[2];
rz(-1.829932) q[3];
sx q[3];
rz(-1.0511755) q[3];
sx q[3];
rz(-1.8478307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0084429) q[0];
sx q[0];
rz(-2.382998) q[0];
sx q[0];
rz(-1.5822509) q[0];
rz(-0.82296952) q[1];
sx q[1];
rz(-2.6428849) q[1];
sx q[1];
rz(-1.71436) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2555851) q[0];
sx q[0];
rz(-2.581121) q[0];
sx q[0];
rz(-2.1917402) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64325579) q[2];
sx q[2];
rz(-1.7417041) q[2];
sx q[2];
rz(-2.3115445) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.848915) q[1];
sx q[1];
rz(-2.1452322) q[1];
sx q[1];
rz(0.019197063) q[1];
rz(-pi) q[2];
rz(0.41540878) q[3];
sx q[3];
rz(-1.6744288) q[3];
sx q[3];
rz(1.8745916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3048749) q[2];
sx q[2];
rz(-2.0676401) q[2];
sx q[2];
rz(2.7433336) q[2];
rz(-0.38175976) q[3];
sx q[3];
rz(-1.7002117) q[3];
sx q[3];
rz(-2.1329342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1036103) q[0];
sx q[0];
rz(-1.4382265) q[0];
sx q[0];
rz(-1.2305228) q[0];
rz(-1.1331753) q[1];
sx q[1];
rz(-0.76825348) q[1];
sx q[1];
rz(2.5450328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.848551) q[0];
sx q[0];
rz(-2.8900091) q[0];
sx q[0];
rz(0.037034794) q[0];
rz(-pi) q[1];
rz(1.0870129) q[2];
sx q[2];
rz(-1.7426335) q[2];
sx q[2];
rz(0.89110723) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0016252) q[1];
sx q[1];
rz(-2.4974065) q[1];
sx q[1];
rz(-2.0124547) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9539709) q[3];
sx q[3];
rz(-2.1541693) q[3];
sx q[3];
rz(-2.4125828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29264984) q[2];
sx q[2];
rz(-2.5299447) q[2];
sx q[2];
rz(1.4673648) q[2];
rz(-0.26662207) q[3];
sx q[3];
rz(-0.98054236) q[3];
sx q[3];
rz(-2.9149741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4416696) q[0];
sx q[0];
rz(-3.141098) q[0];
sx q[0];
rz(1.1215522) q[0];
rz(-2.9660405) q[1];
sx q[1];
rz(-2.4708864) q[1];
sx q[1];
rz(2.6480754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4283912) q[0];
sx q[0];
rz(-1.512171) q[0];
sx q[0];
rz(-1.6423663) q[0];
rz(1.1404806) q[2];
sx q[2];
rz(-0.80289932) q[2];
sx q[2];
rz(-2.0446628) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3596489) q[1];
sx q[1];
rz(-1.8033866) q[1];
sx q[1];
rz(2.1185758) q[1];
x q[2];
rz(3.054627) q[3];
sx q[3];
rz(-2.0425597) q[3];
sx q[3];
rz(1.903002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4404122) q[2];
sx q[2];
rz(-0.36786011) q[2];
sx q[2];
rz(-0.58604798) q[2];
rz(1.4513357) q[3];
sx q[3];
rz(-1.3380932) q[3];
sx q[3];
rz(-0.51226789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.119809) q[0];
sx q[0];
rz(-1.021294) q[0];
sx q[0];
rz(-2.9721066) q[0];
rz(2.4985864) q[1];
sx q[1];
rz(-2.2702859) q[1];
sx q[1];
rz(-0.51173425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6089451) q[0];
sx q[0];
rz(-1.4911312) q[0];
sx q[0];
rz(-0.10124036) q[0];
rz(-2.8795661) q[2];
sx q[2];
rz(-1.507477) q[2];
sx q[2];
rz(-0.89010191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6890672) q[1];
sx q[1];
rz(-2.2241909) q[1];
sx q[1];
rz(-2.9513068) q[1];
rz(-pi) q[2];
rz(0.2318895) q[3];
sx q[3];
rz(-2.0733193) q[3];
sx q[3];
rz(-2.1351867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6181651) q[2];
sx q[2];
rz(-0.7343505) q[2];
sx q[2];
rz(-3.1132474) q[2];
rz(2.667115) q[3];
sx q[3];
rz(-0.35895434) q[3];
sx q[3];
rz(-2.1019782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9356553) q[0];
sx q[0];
rz(-1.7556998) q[0];
sx q[0];
rz(2.6082267) q[0];
rz(-0.43279303) q[1];
sx q[1];
rz(-1.5569893) q[1];
sx q[1];
rz(2.5098262) q[1];
rz(-2.9262056) q[2];
sx q[2];
rz(-1.4949847) q[2];
sx q[2];
rz(1.1204213) q[2];
rz(-3.0632426) q[3];
sx q[3];
rz(-1.363277) q[3];
sx q[3];
rz(-1.1860397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
