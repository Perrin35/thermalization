OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8074789) q[0];
sx q[0];
rz(-0.69598323) q[0];
sx q[0];
rz(2.2544883) q[0];
rz(2.4298985) q[1];
sx q[1];
rz(-2.0790172) q[1];
sx q[1];
rz(0.33023155) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45107156) q[0];
sx q[0];
rz(-2.1981648) q[0];
sx q[0];
rz(2.5894707) q[0];
x q[1];
rz(2.5862972) q[2];
sx q[2];
rz(-1.2568297) q[2];
sx q[2];
rz(1.2295251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.020259) q[1];
sx q[1];
rz(-1.4687045) q[1];
sx q[1];
rz(2.6481204) q[1];
rz(0.078465538) q[3];
sx q[3];
rz(-2.0002504) q[3];
sx q[3];
rz(-1.0965018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3868788) q[2];
sx q[2];
rz(-0.42749307) q[2];
sx q[2];
rz(-2.1145561) q[2];
rz(-0.88296452) q[3];
sx q[3];
rz(-1.3436147) q[3];
sx q[3];
rz(-2.5880255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7850194) q[0];
sx q[0];
rz(-3.0632601) q[0];
sx q[0];
rz(0.26764348) q[0];
rz(-1.74125) q[1];
sx q[1];
rz(-0.61336556) q[1];
sx q[1];
rz(2.5852481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436377) q[0];
sx q[0];
rz(-2.58444) q[0];
sx q[0];
rz(-2.055026) q[0];
rz(0.70425561) q[2];
sx q[2];
rz(-1.8374763) q[2];
sx q[2];
rz(1.2885338) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9832335) q[1];
sx q[1];
rz(-0.69770798) q[1];
sx q[1];
rz(-2.1844728) q[1];
rz(-pi) q[2];
rz(-0.90640599) q[3];
sx q[3];
rz(-1.0944546) q[3];
sx q[3];
rz(-1.9204467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9050425) q[2];
sx q[2];
rz(-1.6908129) q[2];
sx q[2];
rz(1.8897918) q[2];
rz(1.8768138) q[3];
sx q[3];
rz(-0.71139657) q[3];
sx q[3];
rz(2.0701764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65997893) q[0];
sx q[0];
rz(-0.55584207) q[0];
sx q[0];
rz(-0.83928338) q[0];
rz(1.6370157) q[1];
sx q[1];
rz(-1.7103651) q[1];
sx q[1];
rz(-1.379871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312442) q[0];
sx q[0];
rz(-1.6524685) q[0];
sx q[0];
rz(-1.4468589) q[0];
rz(-1.0255541) q[2];
sx q[2];
rz(-1.386706) q[2];
sx q[2];
rz(1.4410401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31095895) q[1];
sx q[1];
rz(-2.0409901) q[1];
sx q[1];
rz(-0.36010919) q[1];
x q[2];
rz(-0.57137892) q[3];
sx q[3];
rz(-1.0927754) q[3];
sx q[3];
rz(-0.92438053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69924361) q[2];
sx q[2];
rz(-1.9244497) q[2];
sx q[2];
rz(2.4231363) q[2];
rz(-0.62026223) q[3];
sx q[3];
rz(-2.2060427) q[3];
sx q[3];
rz(2.3292495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45850596) q[0];
sx q[0];
rz(-1.9106671) q[0];
sx q[0];
rz(1.7339535) q[0];
rz(0.55157026) q[1];
sx q[1];
rz(-0.11081129) q[1];
sx q[1];
rz(1.3161906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39945093) q[0];
sx q[0];
rz(-2.4654498) q[0];
sx q[0];
rz(2.0318982) q[0];
rz(-pi) q[1];
rz(3.1291601) q[2];
sx q[2];
rz(-1.8502478) q[2];
sx q[2];
rz(-0.81571992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0312457) q[1];
sx q[1];
rz(-1.7025885) q[1];
sx q[1];
rz(2.736997) q[1];
rz(-pi) q[2];
rz(0.67263453) q[3];
sx q[3];
rz(-0.97390538) q[3];
sx q[3];
rz(1.2850645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59336415) q[2];
sx q[2];
rz(-2.7061988) q[2];
sx q[2];
rz(1.0080053) q[2];
rz(-2.8625782) q[3];
sx q[3];
rz(-2.3271826) q[3];
sx q[3];
rz(-2.6830955) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27610436) q[0];
sx q[0];
rz(-1.8115598) q[0];
sx q[0];
rz(2.8651067) q[0];
rz(-2.8522988) q[1];
sx q[1];
rz(-1.9464867) q[1];
sx q[1];
rz(2.8791265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55928265) q[0];
sx q[0];
rz(-1.773343) q[0];
sx q[0];
rz(1.7515394) q[0];
x q[1];
rz(0.010774537) q[2];
sx q[2];
rz(-0.97948217) q[2];
sx q[2];
rz(2.7464405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.561132) q[1];
sx q[1];
rz(-1.2314241) q[1];
sx q[1];
rz(2.5575693) q[1];
rz(-pi) q[2];
rz(0.57837242) q[3];
sx q[3];
rz(-1.4422073) q[3];
sx q[3];
rz(-1.8687539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3929954) q[2];
sx q[2];
rz(-0.59610072) q[2];
sx q[2];
rz(-1.7945012) q[2];
rz(-0.10284452) q[3];
sx q[3];
rz(-1.2343854) q[3];
sx q[3];
rz(1.5170521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0031925072) q[0];
sx q[0];
rz(-1.6380558) q[0];
sx q[0];
rz(-3.0806155) q[0];
rz(-1.2334476) q[1];
sx q[1];
rz(-1.7465218) q[1];
sx q[1];
rz(1.5256418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87729851) q[0];
sx q[0];
rz(-1.0420615) q[0];
sx q[0];
rz(1.4905457) q[0];
x q[1];
rz(2.3198177) q[2];
sx q[2];
rz(-1.5944325) q[2];
sx q[2];
rz(1.826959) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96887302) q[1];
sx q[1];
rz(-2.8668826) q[1];
sx q[1];
rz(1.151888) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7876225) q[3];
sx q[3];
rz(-0.66321709) q[3];
sx q[3];
rz(1.4233703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6700651) q[2];
sx q[2];
rz(-1.7933041) q[2];
sx q[2];
rz(-1.8111551) q[2];
rz(1.1250251) q[3];
sx q[3];
rz(-0.87441134) q[3];
sx q[3];
rz(-3.103638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1838945) q[0];
sx q[0];
rz(-1.4426458) q[0];
sx q[0];
rz(-0.25699064) q[0];
rz(-2.7067302) q[1];
sx q[1];
rz(-0.75006524) q[1];
sx q[1];
rz(2.238067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5241476) q[0];
sx q[0];
rz(-3.0492231) q[0];
sx q[0];
rz(1.525982) q[0];
rz(2.3955879) q[2];
sx q[2];
rz(-1.6705319) q[2];
sx q[2];
rz(-1.4745281) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4005247) q[1];
sx q[1];
rz(-1.5092821) q[1];
sx q[1];
rz(-1.5232624) q[1];
rz(-pi) q[2];
rz(-1.0205471) q[3];
sx q[3];
rz(-1.9686832) q[3];
sx q[3];
rz(-0.83845058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8187108) q[2];
sx q[2];
rz(-1.5078397) q[2];
sx q[2];
rz(1.4212849) q[2];
rz(0.70982248) q[3];
sx q[3];
rz(-2.0028508) q[3];
sx q[3];
rz(0.53409725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1555971) q[0];
sx q[0];
rz(-0.43102145) q[0];
sx q[0];
rz(2.0565597) q[0];
rz(0.64708465) q[1];
sx q[1];
rz(-1.9661463) q[1];
sx q[1];
rz(0.064402493) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10391271) q[0];
sx q[0];
rz(-1.2453916) q[0];
sx q[0];
rz(-0.66489451) q[0];
x q[1];
rz(2.2930458) q[2];
sx q[2];
rz(-2.1050958) q[2];
sx q[2];
rz(-1.2017182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29147061) q[1];
sx q[1];
rz(-2.2434592) q[1];
sx q[1];
rz(2.5419278) q[1];
rz(-pi) q[2];
rz(-1.2149947) q[3];
sx q[3];
rz(-1.8321773) q[3];
sx q[3];
rz(-1.3047117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25962466) q[2];
sx q[2];
rz(-0.3825863) q[2];
sx q[2];
rz(-0.043370334) q[2];
rz(-0.67443332) q[3];
sx q[3];
rz(-1.3573656) q[3];
sx q[3];
rz(1.1446713) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944273) q[0];
sx q[0];
rz(-2.5612216) q[0];
sx q[0];
rz(2.8149783) q[0];
rz(-1.0940374) q[1];
sx q[1];
rz(-0.84427619) q[1];
sx q[1];
rz(-0.48386595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1492067) q[0];
sx q[0];
rz(-1.0629553) q[0];
sx q[0];
rz(-0.81224982) q[0];
x q[1];
rz(-2.7174453) q[2];
sx q[2];
rz(-1.3636936) q[2];
sx q[2];
rz(-2.066156) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6623508) q[1];
sx q[1];
rz(-1.3562981) q[1];
sx q[1];
rz(0.7995601) q[1];
x q[2];
rz(-0.30535474) q[3];
sx q[3];
rz(-0.66514665) q[3];
sx q[3];
rz(-2.0572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4491552) q[2];
sx q[2];
rz(-1.038895) q[2];
sx q[2];
rz(0.078744002) q[2];
rz(-2.3732843) q[3];
sx q[3];
rz(-1.8470304) q[3];
sx q[3];
rz(-2.9505742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67097265) q[0];
sx q[0];
rz(-0.60523954) q[0];
sx q[0];
rz(-3.0249) q[0];
rz(2.281588) q[1];
sx q[1];
rz(-1.8000894) q[1];
sx q[1];
rz(2.2681627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1310932) q[0];
sx q[0];
rz(-2.2981055) q[0];
sx q[0];
rz(2.8167679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.248718) q[2];
sx q[2];
rz(-1.0828185) q[2];
sx q[2];
rz(-0.26103324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22120093) q[1];
sx q[1];
rz(-1.4468235) q[1];
sx q[1];
rz(2.6392379) q[1];
rz(-pi) q[2];
rz(0.61140538) q[3];
sx q[3];
rz(-1.5152001) q[3];
sx q[3];
rz(1.1691332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42150911) q[2];
sx q[2];
rz(-2.8953711) q[2];
sx q[2];
rz(-0.949238) q[2];
rz(-1.2468437) q[3];
sx q[3];
rz(-1.6722164) q[3];
sx q[3];
rz(-2.5505572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.128189) q[0];
sx q[0];
rz(-1.0008151) q[0];
sx q[0];
rz(-1.6304954) q[0];
rz(0.36698256) q[1];
sx q[1];
rz(-1.3584247) q[1];
sx q[1];
rz(-0.57986837) q[1];
rz(-1.0147711) q[2];
sx q[2];
rz(-0.53511878) q[2];
sx q[2];
rz(-2.6352885) q[2];
rz(1.9909158) q[3];
sx q[3];
rz(-0.92460604) q[3];
sx q[3];
rz(-2.9210319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
