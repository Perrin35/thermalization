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
rz(-1.3462525) q[0];
sx q[0];
rz(-1.4486382) q[0];
sx q[0];
rz(1.1516655) q[0];
rz(-2.0517853) q[1];
sx q[1];
rz(-1.6648219) q[1];
sx q[1];
rz(2.9906315) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.359084) q[0];
sx q[0];
rz(-1.5182966) q[0];
sx q[0];
rz(-1.8655325) q[0];
rz(-pi) q[1];
rz(1.6577166) q[2];
sx q[2];
rz(-1.8833232) q[2];
sx q[2];
rz(-0.0258044) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0651567) q[1];
sx q[1];
rz(-1.5671092) q[1];
sx q[1];
rz(-1.5855968) q[1];
rz(-3.0449115) q[3];
sx q[3];
rz(-1.5301609) q[3];
sx q[3];
rz(-1.7698897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0207409) q[2];
sx q[2];
rz(-0.023422478) q[2];
sx q[2];
rz(0.49129301) q[2];
rz(1.322572) q[3];
sx q[3];
rz(-1.6072075) q[3];
sx q[3];
rz(-1.8137431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78111929) q[0];
sx q[0];
rz(-2.5871215) q[0];
sx q[0];
rz(2.4419899) q[0];
rz(1.5910925) q[1];
sx q[1];
rz(-0.49483776) q[1];
sx q[1];
rz(2.9717305) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2401597) q[0];
sx q[0];
rz(-0.086726464) q[0];
sx q[0];
rz(-0.46115748) q[0];
x q[1];
rz(1.705395) q[2];
sx q[2];
rz(-1.8730117) q[2];
sx q[2];
rz(-2.5562191) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7034727) q[1];
sx q[1];
rz(-1.2143163) q[1];
sx q[1];
rz(-2.0894024) q[1];
rz(-0.079143123) q[3];
sx q[3];
rz(-2.0266266) q[3];
sx q[3];
rz(-1.7439708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26213172) q[2];
sx q[2];
rz(-2.7792271) q[2];
sx q[2];
rz(1.4169089) q[2];
rz(2.0587685) q[3];
sx q[3];
rz(-2.2542451) q[3];
sx q[3];
rz(2.3143903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.5812434) q[0];
sx q[0];
rz(-0.92009783) q[0];
sx q[0];
rz(-1.5823407) q[0];
rz(2.4371367) q[1];
sx q[1];
rz(-1.1471986) q[1];
sx q[1];
rz(-0.53680435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052633135) q[0];
sx q[0];
rz(-0.68872243) q[0];
sx q[0];
rz(0.072865085) q[0];
rz(-pi) q[1];
rz(0.090254668) q[2];
sx q[2];
rz(-1.5343416) q[2];
sx q[2];
rz(2.0608471) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0672863) q[1];
sx q[1];
rz(-2.9565449) q[1];
sx q[1];
rz(-2.5230899) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87450366) q[3];
sx q[3];
rz(-1.7603959) q[3];
sx q[3];
rz(2.7087871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.392936) q[2];
sx q[2];
rz(-1.8042678) q[2];
sx q[2];
rz(2.362747) q[2];
rz(-3.0574162) q[3];
sx q[3];
rz(-0.86936969) q[3];
sx q[3];
rz(1.7387773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8970784) q[0];
sx q[0];
rz(-0.92105138) q[0];
sx q[0];
rz(-0.47065863) q[0];
rz(-1.4587519) q[1];
sx q[1];
rz(-0.0048480821) q[1];
sx q[1];
rz(0.86476129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5697073) q[0];
sx q[0];
rz(-1.0568468) q[0];
sx q[0];
rz(-1.6374915) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7030348) q[2];
sx q[2];
rz(-1.3515499) q[2];
sx q[2];
rz(-1.8869149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75820747) q[1];
sx q[1];
rz(-1.5050737) q[1];
sx q[1];
rz(-0.051082146) q[1];
rz(-pi) q[2];
rz(-0.49089499) q[3];
sx q[3];
rz(-0.59661907) q[3];
sx q[3];
rz(-0.20363775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.027815759) q[2];
sx q[2];
rz(-2.2769589) q[2];
sx q[2];
rz(1.1569542) q[2];
rz(0.35465869) q[3];
sx q[3];
rz(-1.9054474) q[3];
sx q[3];
rz(-3.0310596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81195867) q[0];
sx q[0];
rz(-0.93557731) q[0];
sx q[0];
rz(-2.6079566) q[0];
rz(-1.927902) q[1];
sx q[1];
rz(-0.023093725) q[1];
sx q[1];
rz(1.9603221) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.597724) q[0];
sx q[0];
rz(-1.7337611) q[0];
sx q[0];
rz(-1.5242432) q[0];
rz(-2.6152739) q[2];
sx q[2];
rz(-0.72229702) q[2];
sx q[2];
rz(-0.96529759) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3141796) q[1];
sx q[1];
rz(-1.7269969) q[1];
sx q[1];
rz(-0.91663702) q[1];
rz(-1.817068) q[3];
sx q[3];
rz(-0.51760537) q[3];
sx q[3];
rz(0.31806614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10826762) q[2];
sx q[2];
rz(-2.6259618) q[2];
sx q[2];
rz(0.92833129) q[2];
rz(-0.27836529) q[3];
sx q[3];
rz(-1.3185578) q[3];
sx q[3];
rz(-3.0729455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6081029) q[0];
sx q[0];
rz(-0.52284306) q[0];
sx q[0];
rz(0.19304092) q[0];
rz(0.67735425) q[1];
sx q[1];
rz(-3.1132128) q[1];
sx q[1];
rz(1.631564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0583147) q[0];
sx q[0];
rz(-1.4299576) q[0];
sx q[0];
rz(1.550564) q[0];
x q[1];
rz(-1.6608428) q[2];
sx q[2];
rz(-2.6777732) q[2];
sx q[2];
rz(0.12741379) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2605156) q[1];
sx q[1];
rz(-2.279638) q[1];
sx q[1];
rz(-2.8167732) q[1];
rz(2.4149309) q[3];
sx q[3];
rz(-1.1805941) q[3];
sx q[3];
rz(1.778804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3626927) q[2];
sx q[2];
rz(-0.71234667) q[2];
sx q[2];
rz(-2.1544797) q[2];
rz(2.0569862) q[3];
sx q[3];
rz(-0.54607138) q[3];
sx q[3];
rz(2.5206821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5622332) q[0];
sx q[0];
rz(-2.1491304) q[0];
sx q[0];
rz(1.5935422) q[0];
rz(-0.69001895) q[1];
sx q[1];
rz(-0.023612173) q[1];
sx q[1];
rz(1.7049559) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6520335) q[0];
sx q[0];
rz(-1.0512182) q[0];
sx q[0];
rz(1.5205617) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1146554) q[2];
sx q[2];
rz(-0.68490324) q[2];
sx q[2];
rz(1.2482582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0322528) q[1];
sx q[1];
rz(-0.54530646) q[1];
sx q[1];
rz(1.6857573) q[1];
rz(-1.8687181) q[3];
sx q[3];
rz(-1.1323338) q[3];
sx q[3];
rz(2.7496453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7684266) q[2];
sx q[2];
rz(-1.2370141) q[2];
sx q[2];
rz(-0.58488673) q[2];
rz(-1.54555) q[3];
sx q[3];
rz(-1.7187748) q[3];
sx q[3];
rz(-0.91442951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5698513) q[0];
sx q[0];
rz(-0.91747704) q[0];
sx q[0];
rz(-1.566994) q[0];
rz(2.6245116) q[1];
sx q[1];
rz(-2.2061429) q[1];
sx q[1];
rz(1.2665952) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168473) q[0];
sx q[0];
rz(-0.006082919) q[0];
sx q[0];
rz(-2.541126) q[0];
rz(-pi) q[1];
rz(-0.53740737) q[2];
sx q[2];
rz(-2.2610342) q[2];
sx q[2];
rz(1.7459401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.32990202) q[1];
sx q[1];
rz(-1.2214097) q[1];
sx q[1];
rz(2.1495649) q[1];
rz(1.6932159) q[3];
sx q[3];
rz(-1.4493692) q[3];
sx q[3];
rz(1.1894912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.37303415) q[2];
sx q[2];
rz(-0.80058241) q[2];
sx q[2];
rz(-1.2726834) q[2];
rz(-0.4396762) q[3];
sx q[3];
rz(-1.1594073) q[3];
sx q[3];
rz(-3.0769949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8398297) q[0];
sx q[0];
rz(-3.1253212) q[0];
sx q[0];
rz(-2.8453258) q[0];
rz(0.06427327) q[1];
sx q[1];
rz(-1.4646894) q[1];
sx q[1];
rz(1.6305264) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8730174) q[0];
sx q[0];
rz(-0.35726705) q[0];
sx q[0];
rz(-1.1910458) q[0];
rz(-pi) q[1];
rz(-2.4687541) q[2];
sx q[2];
rz(-1.4110067) q[2];
sx q[2];
rz(0.21110134) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4402995) q[1];
sx q[1];
rz(-1.2475752) q[1];
sx q[1];
rz(-3.0053291) q[1];
rz(-3.1301163) q[3];
sx q[3];
rz(-1.6919129) q[3];
sx q[3];
rz(-2.3055784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5955547) q[2];
sx q[2];
rz(-1.7067355) q[2];
sx q[2];
rz(-1.1240553) q[2];
rz(-1.3042287) q[3];
sx q[3];
rz(-0.29574695) q[3];
sx q[3];
rz(3.0834294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703363) q[0];
sx q[0];
rz(-0.79126343) q[0];
sx q[0];
rz(1.9747718) q[0];
rz(1.600949) q[1];
sx q[1];
rz(-0.29778844) q[1];
sx q[1];
rz(1.8150394) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013915283) q[0];
sx q[0];
rz(-1.3260889) q[0];
sx q[0];
rz(-2.4726333) q[0];
x q[1];
rz(2.7339513) q[2];
sx q[2];
rz(-1.2349814) q[2];
sx q[2];
rz(-1.6419572) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5759144) q[1];
sx q[1];
rz(-0.88972461) q[1];
sx q[1];
rz(-1.2965864) q[1];
rz(-pi) q[2];
rz(1.79779) q[3];
sx q[3];
rz(-0.41068893) q[3];
sx q[3];
rz(2.3410377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2549609) q[2];
sx q[2];
rz(-0.16356629) q[2];
sx q[2];
rz(-1.6939885) q[2];
rz(3.0149031) q[3];
sx q[3];
rz(-1.6619752) q[3];
sx q[3];
rz(2.0255585) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806725) q[0];
sx q[0];
rz(-1.3074449) q[0];
sx q[0];
rz(1.773651) q[0];
rz(-1.5927636) q[1];
sx q[1];
rz(-0.82850414) q[1];
sx q[1];
rz(-2.9864476) q[1];
rz(2.6947179) q[2];
sx q[2];
rz(-1.7394273) q[2];
sx q[2];
rz(1.5706825) q[2];
rz(-0.98649377) q[3];
sx q[3];
rz(-1.7135847) q[3];
sx q[3];
rz(1.8691487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
