OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4898981) q[0];
sx q[0];
rz(-2.2597921) q[0];
sx q[0];
rz(0.14468004) q[0];
rz(3.559685) q[1];
sx q[1];
rz(3.2453645) q[1];
sx q[1];
rz(11.054463) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3016925) q[0];
sx q[0];
rz(-1.6891251) q[0];
sx q[0];
rz(-2.0685643) q[0];
x q[1];
rz(-2.4080963) q[2];
sx q[2];
rz(-1.6111177) q[2];
sx q[2];
rz(0.87903838) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6408893) q[1];
sx q[1];
rz(-1.2391866) q[1];
sx q[1];
rz(-1.160781) q[1];
rz(1.9910664) q[3];
sx q[3];
rz(-2.0448677) q[3];
sx q[3];
rz(-1.898511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8251553) q[2];
sx q[2];
rz(-1.2707767) q[2];
sx q[2];
rz(2.4862508) q[2];
rz(-1.7573028) q[3];
sx q[3];
rz(-2.1770848) q[3];
sx q[3];
rz(0.82096076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69219387) q[0];
sx q[0];
rz(-1.7658424) q[0];
sx q[0];
rz(-0.50814116) q[0];
rz(-1.3805768) q[1];
sx q[1];
rz(-2.5602129) q[1];
sx q[1];
rz(1.2095721) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.094645) q[0];
sx q[0];
rz(-2.1535465) q[0];
sx q[0];
rz(0.21999448) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0879966) q[2];
sx q[2];
rz(-1.5884678) q[2];
sx q[2];
rz(-2.0272209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3833419) q[1];
sx q[1];
rz(-2.0313325) q[1];
sx q[1];
rz(1.7136445) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1749635) q[3];
sx q[3];
rz(-1.2467017) q[3];
sx q[3];
rz(-1.54502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2035344) q[2];
sx q[2];
rz(-1.7455696) q[2];
sx q[2];
rz(2.5246942) q[2];
rz(-1.5361702) q[3];
sx q[3];
rz(-2.0974396) q[3];
sx q[3];
rz(0.49083403) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598375) q[0];
sx q[0];
rz(-1.6438537) q[0];
sx q[0];
rz(2.4165261) q[0];
rz(1.8720576) q[1];
sx q[1];
rz(-0.67765647) q[1];
sx q[1];
rz(-0.95917541) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0051668) q[0];
sx q[0];
rz(-1.1014897) q[0];
sx q[0];
rz(-2.704436) q[0];
x q[1];
rz(-0.74752918) q[2];
sx q[2];
rz(-0.5420712) q[2];
sx q[2];
rz(-2.4666748) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37466418) q[1];
sx q[1];
rz(-0.21682021) q[1];
sx q[1];
rz(2.1310852) q[1];
x q[2];
rz(-1.3884344) q[3];
sx q[3];
rz(-1.404247) q[3];
sx q[3];
rz(0.20660755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3265257) q[2];
sx q[2];
rz(-1.2667789) q[2];
sx q[2];
rz(-0.25700021) q[2];
rz(-1.0810931) q[3];
sx q[3];
rz(-2.6205781) q[3];
sx q[3];
rz(-1.5787554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0576039) q[0];
sx q[0];
rz(-0.27844089) q[0];
sx q[0];
rz(1.9637015) q[0];
rz(2.9914757) q[1];
sx q[1];
rz(-0.53214407) q[1];
sx q[1];
rz(1.6252801) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0310136) q[0];
sx q[0];
rz(-0.4415126) q[0];
sx q[0];
rz(-2.2792321) q[0];
rz(-pi) q[1];
rz(0.84211911) q[2];
sx q[2];
rz(-1.2884022) q[2];
sx q[2];
rz(0.51674622) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8575688) q[1];
sx q[1];
rz(-2.800436) q[1];
sx q[1];
rz(-0.30156231) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.890804) q[3];
sx q[3];
rz(-0.65902519) q[3];
sx q[3];
rz(3.0626754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2991422) q[2];
sx q[2];
rz(-2.5717042) q[2];
sx q[2];
rz(-0.92310706) q[2];
rz(2.182377) q[3];
sx q[3];
rz(-2.1289289) q[3];
sx q[3];
rz(-0.21624163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8206896) q[0];
sx q[0];
rz(-0.90383363) q[0];
sx q[0];
rz(-0.88291105) q[0];
rz(1.2954905) q[1];
sx q[1];
rz(-0.77406445) q[1];
sx q[1];
rz(-0.49497089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6225016) q[0];
sx q[0];
rz(-1.8880196) q[0];
sx q[0];
rz(-2.836801) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4385543) q[2];
sx q[2];
rz(-1.993317) q[2];
sx q[2];
rz(-0.51681821) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.043083) q[1];
sx q[1];
rz(-1.483023) q[1];
sx q[1];
rz(-0.66941525) q[1];
rz(-pi) q[2];
rz(-0.35843973) q[3];
sx q[3];
rz(-2.0652986) q[3];
sx q[3];
rz(0.93264893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3317269) q[2];
sx q[2];
rz(-3.0220384) q[2];
sx q[2];
rz(-1.2681819) q[2];
rz(-0.082503334) q[3];
sx q[3];
rz(-1.637633) q[3];
sx q[3];
rz(0.65703195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6562011) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(1.50151) q[0];
rz(-2.0791176) q[1];
sx q[1];
rz(-1.1524009) q[1];
sx q[1];
rz(1.1753488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6268183) q[0];
sx q[0];
rz(-1.8747703) q[0];
sx q[0];
rz(1.8016812) q[0];
rz(-pi) q[1];
rz(-2.7926867) q[2];
sx q[2];
rz(-2.4173173) q[2];
sx q[2];
rz(1.5021715) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.7580686) q[1];
sx q[1];
rz(-1.3192156) q[1];
sx q[1];
rz(-1.6918159) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31194056) q[3];
sx q[3];
rz(-0.7850724) q[3];
sx q[3];
rz(-1.5573481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33092734) q[2];
sx q[2];
rz(-1.284548) q[2];
sx q[2];
rz(2.1539099) q[2];
rz(2.2641613) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(-0.73006829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1838609) q[0];
sx q[0];
rz(-1.6022302) q[0];
sx q[0];
rz(-3.0027332) q[0];
rz(1.214437) q[1];
sx q[1];
rz(-1.5077533) q[1];
sx q[1];
rz(2.3249998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2557295) q[0];
sx q[0];
rz(-0.22561377) q[0];
sx q[0];
rz(1.4855235) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8726878) q[2];
sx q[2];
rz(-0.78654002) q[2];
sx q[2];
rz(-1.5812909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6703709) q[1];
sx q[1];
rz(-2.2891217) q[1];
sx q[1];
rz(-0.049777546) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2995629) q[3];
sx q[3];
rz(-0.86787628) q[3];
sx q[3];
rz(-1.7563535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.82178086) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(0.4846586) q[2];
rz(-2.7759077) q[3];
sx q[3];
rz(-1.7771143) q[3];
sx q[3];
rz(1.1588089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1300238) q[0];
sx q[0];
rz(-0.31112177) q[0];
sx q[0];
rz(0.097231641) q[0];
rz(-3.0367127) q[1];
sx q[1];
rz(-0.77969867) q[1];
sx q[1];
rz(-1.3146776) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0785694) q[0];
sx q[0];
rz(-0.0067575909) q[0];
sx q[0];
rz(-1.2122173) q[0];
x q[1];
rz(-1.8123367) q[2];
sx q[2];
rz(-1.7881696) q[2];
sx q[2];
rz(-1.5940983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2995616) q[1];
sx q[1];
rz(-0.52134575) q[1];
sx q[1];
rz(-0.10592769) q[1];
rz(-pi) q[2];
rz(-2.1076074) q[3];
sx q[3];
rz(-1.4958463) q[3];
sx q[3];
rz(2.1225464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.17681992) q[2];
sx q[2];
rz(-1.7440045) q[2];
sx q[2];
rz(-2.9883265) q[2];
rz(-1.2166474) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(-2.5254068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0230873) q[0];
sx q[0];
rz(-3.1361134) q[0];
sx q[0];
rz(1.6368921) q[0];
rz(-0.79832375) q[1];
sx q[1];
rz(-1.5232122) q[1];
sx q[1];
rz(2.8062779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3753877) q[0];
sx q[0];
rz(-1.5076625) q[0];
sx q[0];
rz(1.8324018) q[0];
x q[1];
rz(1.5952571) q[2];
sx q[2];
rz(-2.3210397) q[2];
sx q[2];
rz(0.94708196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0287231) q[1];
sx q[1];
rz(-1.8577033) q[1];
sx q[1];
rz(2.7863281) q[1];
rz(-pi) q[2];
rz(-0.15729372) q[3];
sx q[3];
rz(-0.99911896) q[3];
sx q[3];
rz(-1.2441105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12751427) q[2];
sx q[2];
rz(-1.9731584) q[2];
sx q[2];
rz(-2.7039995) q[2];
rz(-1.0420927) q[3];
sx q[3];
rz(-1.4053248) q[3];
sx q[3];
rz(-2.5991345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1522778) q[0];
sx q[0];
rz(-2.2570026) q[0];
sx q[0];
rz(0.44961318) q[0];
rz(0.97688976) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(-2.6403715) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6527443) q[0];
sx q[0];
rz(-2.7978659) q[0];
sx q[0];
rz(2.1656242) q[0];
rz(2.0222763) q[2];
sx q[2];
rz(-0.46487936) q[2];
sx q[2];
rz(1.7865739) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3851191) q[1];
sx q[1];
rz(-1.8071399) q[1];
sx q[1];
rz(2.0467351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1057749) q[3];
sx q[3];
rz(-1.8793725) q[3];
sx q[3];
rz(3.05911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0235128) q[2];
sx q[2];
rz(-2.2492275) q[2];
sx q[2];
rz(0.21772131) q[2];
rz(1.2348385) q[3];
sx q[3];
rz(-2.4775938) q[3];
sx q[3];
rz(-0.92065221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49391838) q[0];
sx q[0];
rz(-1.7128581) q[0];
sx q[0];
rz(0.38059522) q[0];
rz(2.4299798) q[1];
sx q[1];
rz(-1.2672392) q[1];
sx q[1];
rz(-1.738501) q[1];
rz(-0.099945036) q[2];
sx q[2];
rz(-0.47158416) q[2];
sx q[2];
rz(-0.33000962) q[2];
rz(-1.8071411) q[3];
sx q[3];
rz(-1.1750887) q[3];
sx q[3];
rz(-0.25050161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
