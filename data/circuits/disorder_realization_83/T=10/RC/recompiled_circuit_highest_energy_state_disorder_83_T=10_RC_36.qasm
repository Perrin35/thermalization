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
rz(2.7163765) q[0];
sx q[0];
rz(-1.334231) q[0];
sx q[0];
rz(-2.7612326) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(-0.18667297) q[1];
sx q[1];
rz(-2.2169854) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0565489) q[0];
sx q[0];
rz(-0.88221545) q[0];
sx q[0];
rz(-1.1231281) q[0];
x q[1];
rz(2.9745462) q[2];
sx q[2];
rz(-1.0606597) q[2];
sx q[2];
rz(1.270592) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8421665) q[1];
sx q[1];
rz(-1.5314845) q[1];
sx q[1];
rz(-2.2996603) q[1];
x q[2];
rz(-0.070234039) q[3];
sx q[3];
rz(-0.12824225) q[3];
sx q[3];
rz(-1.8907036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1251462) q[2];
sx q[2];
rz(-2.1852198) q[2];
sx q[2];
rz(2.1535786) q[2];
rz(0.2068578) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(-0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67966953) q[0];
sx q[0];
rz(-1.672687) q[0];
sx q[0];
rz(0.2521421) q[0];
rz(0.02154669) q[1];
sx q[1];
rz(-2.4485059) q[1];
sx q[1];
rz(-2.2861939) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48777521) q[0];
sx q[0];
rz(-3.1145373) q[0];
sx q[0];
rz(1.5144801) q[0];
rz(-pi) q[1];
rz(-2.5163745) q[2];
sx q[2];
rz(-0.65186912) q[2];
sx q[2];
rz(-0.21833459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63355061) q[1];
sx q[1];
rz(-1.3505858) q[1];
sx q[1];
rz(-1.0255662) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9048433) q[3];
sx q[3];
rz(-0.91558394) q[3];
sx q[3];
rz(1.987059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1713193) q[2];
sx q[2];
rz(-0.0042985175) q[2];
sx q[2];
rz(0.26101905) q[2];
rz(-0.039693443) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(-2.5980914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71565851) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(0.69450992) q[0];
rz(-2.014324) q[1];
sx q[1];
rz(-0.85113168) q[1];
sx q[1];
rz(0.99004254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111525) q[0];
sx q[0];
rz(-1.2227204) q[0];
sx q[0];
rz(-1.6578107) q[0];
x q[1];
rz(-1.8672529) q[2];
sx q[2];
rz(-0.60288376) q[2];
sx q[2];
rz(-2.4566513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5824008) q[1];
sx q[1];
rz(-1.2123931) q[1];
sx q[1];
rz(2.6646975) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86750908) q[3];
sx q[3];
rz(-1.1425619) q[3];
sx q[3];
rz(-3.0813062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29828829) q[2];
sx q[2];
rz(-2.1683606) q[2];
sx q[2];
rz(2.6514371) q[2];
rz(0.35689029) q[3];
sx q[3];
rz(-0.39339742) q[3];
sx q[3];
rz(-1.2662158) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056034293) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(-0.84247843) q[0];
rz(0.8017686) q[1];
sx q[1];
rz(-2.8623878) q[1];
sx q[1];
rz(-2.3559949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7467037) q[0];
sx q[0];
rz(-0.75468844) q[0];
sx q[0];
rz(-0.37209038) q[0];
rz(1.4198205) q[2];
sx q[2];
rz(-2.2822126) q[2];
sx q[2];
rz(-2.6670418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.848104) q[1];
sx q[1];
rz(-1.2624718) q[1];
sx q[1];
rz(2.5749194) q[1];
x q[2];
rz(-0.20468851) q[3];
sx q[3];
rz(-0.83016268) q[3];
sx q[3];
rz(1.1799174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12513146) q[2];
sx q[2];
rz(-0.74840122) q[2];
sx q[2];
rz(2.1330736) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9652902) q[0];
sx q[0];
rz(-0.98133123) q[0];
sx q[0];
rz(-2.0710131) q[0];
rz(-2.3640682) q[1];
sx q[1];
rz(-0.91064015) q[1];
sx q[1];
rz(-1.0106962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1436737) q[0];
sx q[0];
rz(-0.61286139) q[0];
sx q[0];
rz(0.65632239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7456876) q[2];
sx q[2];
rz(-1.5064459) q[2];
sx q[2];
rz(1.627587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7091211) q[1];
sx q[1];
rz(-2.2975031) q[1];
sx q[1];
rz(0.62950397) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4713953) q[3];
sx q[3];
rz(-2.1306921) q[3];
sx q[3];
rz(0.82370629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30194482) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(-1.5860175) q[2];
rz(1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(-1.7293845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79065901) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(-0.80107981) q[0];
rz(3.0084897) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(0.4020234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19045705) q[0];
sx q[0];
rz(-1.2371089) q[0];
sx q[0];
rz(2.3526089) q[0];
rz(-1.5779675) q[2];
sx q[2];
rz(-2.8507887) q[2];
sx q[2];
rz(0.70722843) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.0084486246) q[1];
sx q[1];
rz(-2.5767972) q[1];
sx q[1];
rz(-1.0393591) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9803863) q[3];
sx q[3];
rz(-2.4771002) q[3];
sx q[3];
rz(-0.67918577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5037527) q[2];
sx q[2];
rz(-1.4053586) q[2];
sx q[2];
rz(-2.3657738) q[2];
rz(-1.6860298) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(2.1337401) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044947226) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(2.4626379) q[0];
rz(1.8567765) q[1];
sx q[1];
rz(-2.4285451) q[1];
sx q[1];
rz(1.3791893) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75684568) q[0];
sx q[0];
rz(-2.3492536) q[0];
sx q[0];
rz(-1.8154665) q[0];
x q[1];
rz(-1.2913338) q[2];
sx q[2];
rz(-1.6510291) q[2];
sx q[2];
rz(0.81278518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1393362) q[1];
sx q[1];
rz(-0.69076194) q[1];
sx q[1];
rz(0.73378566) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3016765) q[3];
sx q[3];
rz(-2.9342071) q[3];
sx q[3];
rz(-3.0548422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3369559) q[2];
sx q[2];
rz(-0.75504428) q[2];
sx q[2];
rz(1.2389368) q[2];
rz(-1.3321446) q[3];
sx q[3];
rz(-0.76056162) q[3];
sx q[3];
rz(-0.24470394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.1661538) q[0];
sx q[0];
rz(-0.13846692) q[0];
rz(-2.9578517) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(0.26652452) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4999102) q[0];
sx q[0];
rz(-1.2453834) q[0];
sx q[0];
rz(2.9703559) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1600003) q[2];
sx q[2];
rz(-2.5345384) q[2];
sx q[2];
rz(-2.0701054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9723818) q[1];
sx q[1];
rz(-1.5145497) q[1];
sx q[1];
rz(-1.1509658) q[1];
rz(-pi) q[2];
rz(-3.1317461) q[3];
sx q[3];
rz(-2.4146663) q[3];
sx q[3];
rz(-1.3971412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0515685) q[2];
sx q[2];
rz(-1.2568306) q[2];
sx q[2];
rz(0.23452342) q[2];
rz(-1.4759493) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(-0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.804857) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(2.3418703) q[0];
rz(-0.58894482) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(-0.2074997) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070415592) q[0];
sx q[0];
rz(-1.9635597) q[0];
sx q[0];
rz(-1.7123187) q[0];
rz(-pi) q[1];
rz(0.3090119) q[2];
sx q[2];
rz(-2.2591619) q[2];
sx q[2];
rz(2.7613044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60002414) q[1];
sx q[1];
rz(-2.3316521) q[1];
sx q[1];
rz(-1.0066973) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5141684) q[3];
sx q[3];
rz(-2.671173) q[3];
sx q[3];
rz(-2.3868167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3346682) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(-0.36726382) q[2];
rz(1.7043097) q[3];
sx q[3];
rz(-1.1742914) q[3];
sx q[3];
rz(1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87840286) q[0];
sx q[0];
rz(-1.0028361) q[0];
sx q[0];
rz(-0.81018418) q[0];
rz(1.1148249) q[1];
sx q[1];
rz(-0.38433847) q[1];
sx q[1];
rz(-2.5883163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081757717) q[0];
sx q[0];
rz(-2.5011269) q[0];
sx q[0];
rz(-1.6637131) q[0];
rz(-pi) q[1];
rz(0.61364321) q[2];
sx q[2];
rz(-0.89362234) q[2];
sx q[2];
rz(2.8751862) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3685963) q[1];
sx q[1];
rz(-2.8008411) q[1];
sx q[1];
rz(0.98253886) q[1];
rz(-pi) q[2];
rz(-0.35240473) q[3];
sx q[3];
rz(-0.81768546) q[3];
sx q[3];
rz(-1.529983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0193923) q[2];
sx q[2];
rz(-2.4330752) q[2];
sx q[2];
rz(0.23966399) q[2];
rz(-1.3278809) q[3];
sx q[3];
rz(-1.0646822) q[3];
sx q[3];
rz(0.46752587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620621) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(-2.486034) q[1];
sx q[1];
rz(-1.1759023) q[1];
sx q[1];
rz(-2.1942153) q[1];
rz(-2.4782933) q[2];
sx q[2];
rz(-0.51021432) q[2];
sx q[2];
rz(-2.9205657) q[2];
rz(-2.5611193) q[3];
sx q[3];
rz(-1.3571285) q[3];
sx q[3];
rz(0.50783689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
