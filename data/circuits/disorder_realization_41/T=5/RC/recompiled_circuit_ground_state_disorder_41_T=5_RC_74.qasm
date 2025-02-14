OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5858894) q[0];
sx q[0];
rz(-1.1893505) q[0];
sx q[0];
rz(1.9777634) q[0];
rz(2.3948506) q[1];
sx q[1];
rz(-2.876694) q[1];
sx q[1];
rz(-0.17118153) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88306016) q[0];
sx q[0];
rz(-0.2685606) q[0];
sx q[0];
rz(2.1826151) q[0];
x q[1];
rz(-2.2581882) q[2];
sx q[2];
rz(-2.2210178) q[2];
sx q[2];
rz(-1.6028849) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.491558) q[1];
sx q[1];
rz(-1.2753873) q[1];
sx q[1];
rz(-0.00062382767) q[1];
rz(-0.31193093) q[3];
sx q[3];
rz(-0.017695713) q[3];
sx q[3];
rz(-1.9809256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.050194) q[2];
sx q[2];
rz(-1.5768496) q[2];
sx q[2];
rz(-1.0696627) q[2];
rz(-1.4403249) q[3];
sx q[3];
rz(-2.1170719) q[3];
sx q[3];
rz(1.1721771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1759724) q[0];
sx q[0];
rz(-1.4485285) q[0];
sx q[0];
rz(2.5536221) q[0];
rz(-1.2929471) q[1];
sx q[1];
rz(-2.1473532) q[1];
sx q[1];
rz(2.9867244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3414587) q[0];
sx q[0];
rz(-1.4636369) q[0];
sx q[0];
rz(-0.5002621) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6011818) q[2];
sx q[2];
rz(-0.27920461) q[2];
sx q[2];
rz(1.0591194) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.14531222) q[1];
sx q[1];
rz(-0.59050817) q[1];
sx q[1];
rz(2.202321) q[1];
rz(-1.7072466) q[3];
sx q[3];
rz(-1.3605355) q[3];
sx q[3];
rz(-0.36255676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64314848) q[2];
sx q[2];
rz(-0.73068205) q[2];
sx q[2];
rz(-3.0291962) q[2];
rz(0.82320881) q[3];
sx q[3];
rz(-1.221849) q[3];
sx q[3];
rz(-2.6196151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68558973) q[0];
sx q[0];
rz(-2.0528448) q[0];
sx q[0];
rz(0.96813694) q[0];
rz(2.0227506) q[1];
sx q[1];
rz(-1.6066931) q[1];
sx q[1];
rz(1.3333837) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0007933) q[0];
sx q[0];
rz(-2.5171772) q[0];
sx q[0];
rz(-1.8483759) q[0];
rz(2.5081464) q[2];
sx q[2];
rz(-2.434133) q[2];
sx q[2];
rz(-2.5587669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3231877) q[1];
sx q[1];
rz(-1.0099704) q[1];
sx q[1];
rz(-3.0238999) q[1];
rz(-1.2127688) q[3];
sx q[3];
rz(-1.6505401) q[3];
sx q[3];
rz(-0.84313508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94990388) q[2];
sx q[2];
rz(-1.0368985) q[2];
sx q[2];
rz(2.8766768) q[2];
rz(-2.8252937) q[3];
sx q[3];
rz(-2.8079872) q[3];
sx q[3];
rz(1.7501638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3259657) q[0];
sx q[0];
rz(-2.9209324) q[0];
sx q[0];
rz(1.6478446) q[0];
rz(1.4119459) q[1];
sx q[1];
rz(-1.0271881) q[1];
sx q[1];
rz(-2.484201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7058168) q[0];
sx q[0];
rz(-3.0141493) q[0];
sx q[0];
rz(0.61690046) q[0];
rz(-pi) q[1];
rz(-2.3573586) q[2];
sx q[2];
rz(-2.0080697) q[2];
sx q[2];
rz(2.0511049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.58739793) q[1];
sx q[1];
rz(-2.9430591) q[1];
sx q[1];
rz(-1.3819329) q[1];
x q[2];
rz(0.811565) q[3];
sx q[3];
rz(-0.92716427) q[3];
sx q[3];
rz(-1.0652468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21939453) q[2];
sx q[2];
rz(-1.9191091) q[2];
sx q[2];
rz(2.738319) q[2];
rz(1.9379617) q[3];
sx q[3];
rz(-1.6324685) q[3];
sx q[3];
rz(0.50886124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4598684) q[0];
sx q[0];
rz(-0.42400703) q[0];
sx q[0];
rz(-2.5176609) q[0];
rz(-0.722018) q[1];
sx q[1];
rz(-1.7941509) q[1];
sx q[1];
rz(3.0375979) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.390966) q[0];
sx q[0];
rz(-1.4292882) q[0];
sx q[0];
rz(-2.4510989) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9172759) q[2];
sx q[2];
rz(-1.2305546) q[2];
sx q[2];
rz(2.1635596) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6045058) q[1];
sx q[1];
rz(-2.1963122) q[1];
sx q[1];
rz(0.64670678) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4071619) q[3];
sx q[3];
rz(-1.0356257) q[3];
sx q[3];
rz(1.9246235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35088745) q[2];
sx q[2];
rz(-0.71983379) q[2];
sx q[2];
rz(-0.78901115) q[2];
rz(1.2645432) q[3];
sx q[3];
rz(-1.0435373) q[3];
sx q[3];
rz(0.98715442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8378545) q[0];
sx q[0];
rz(-2.4202977) q[0];
sx q[0];
rz(-2.9413132) q[0];
rz(1.6750977) q[1];
sx q[1];
rz(-1.3267582) q[1];
sx q[1];
rz(1.5178348) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0132594) q[0];
sx q[0];
rz(-0.64494123) q[0];
sx q[0];
rz(2.5027184) q[0];
x q[1];
rz(1.8046204) q[2];
sx q[2];
rz(-2.6262198) q[2];
sx q[2];
rz(-2.963954) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38963803) q[1];
sx q[1];
rz(-1.5672958) q[1];
sx q[1];
rz(-2.5506427) q[1];
x q[2];
rz(2.3430941) q[3];
sx q[3];
rz(-1.3991465) q[3];
sx q[3];
rz(-2.5018179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62310654) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(0.34131193) q[2];
rz(-2.2845279) q[3];
sx q[3];
rz(-1.2271481) q[3];
sx q[3];
rz(-2.6763693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.051006) q[0];
sx q[0];
rz(-2.107928) q[0];
sx q[0];
rz(-1.8940014) q[0];
rz(3.0233851) q[1];
sx q[1];
rz(-1.2756313) q[1];
sx q[1];
rz(-0.035621312) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3105594) q[0];
sx q[0];
rz(-0.11971029) q[0];
sx q[0];
rz(1.7605147) q[0];
rz(-2.9526677) q[2];
sx q[2];
rz(-1.9686724) q[2];
sx q[2];
rz(2.973345) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25917398) q[1];
sx q[1];
rz(-0.66773623) q[1];
sx q[1];
rz(0.40294934) q[1];
rz(-pi) q[2];
rz(1.4485967) q[3];
sx q[3];
rz(-2.3482493) q[3];
sx q[3];
rz(0.47652188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9901765) q[2];
sx q[2];
rz(-1.6665062) q[2];
sx q[2];
rz(1.0756294) q[2];
rz(-0.52078024) q[3];
sx q[3];
rz(-2.5213089) q[3];
sx q[3];
rz(-1.5889997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90466475) q[0];
sx q[0];
rz(-1.0228782) q[0];
sx q[0];
rz(-0.9032332) q[0];
rz(1.5029933) q[1];
sx q[1];
rz(-1.4475977) q[1];
sx q[1];
rz(1.2618056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0340817) q[0];
sx q[0];
rz(-0.26824646) q[0];
sx q[0];
rz(-2.111042) q[0];
rz(-pi) q[1];
rz(0.52623279) q[2];
sx q[2];
rz(-1.1723435) q[2];
sx q[2];
rz(-0.72771074) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7664396) q[1];
sx q[1];
rz(-0.10442153) q[1];
sx q[1];
rz(-2.6344643) q[1];
x q[2];
rz(-0.66916211) q[3];
sx q[3];
rz(-0.72297308) q[3];
sx q[3];
rz(-0.97217901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5624076) q[2];
sx q[2];
rz(-2.4895442) q[2];
sx q[2];
rz(-0.11422608) q[2];
rz(-1.3769897) q[3];
sx q[3];
rz(-2.0012794) q[3];
sx q[3];
rz(0.50447869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1219516) q[0];
sx q[0];
rz(-2.1421102) q[0];
sx q[0];
rz(1.9695388) q[0];
rz(1.7578341) q[1];
sx q[1];
rz(-1.1446605) q[1];
sx q[1];
rz(3.0269472) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5929554) q[0];
sx q[0];
rz(-2.1229738) q[0];
sx q[0];
rz(-2.1594285) q[0];
rz(-2.9681712) q[2];
sx q[2];
rz(-0.49931128) q[2];
sx q[2];
rz(-0.87022802) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.35797) q[1];
sx q[1];
rz(-1.8183892) q[1];
sx q[1];
rz(0.072958306) q[1];
x q[2];
rz(-0.10520868) q[3];
sx q[3];
rz(-2.3423839) q[3];
sx q[3];
rz(1.5952283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5961479) q[2];
sx q[2];
rz(-0.20716509) q[2];
sx q[2];
rz(-0.88085112) q[2];
rz(0.43978575) q[3];
sx q[3];
rz(-1.9763016) q[3];
sx q[3];
rz(-2.0161276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3038444) q[0];
sx q[0];
rz(-1.6551908) q[0];
sx q[0];
rz(0.098175511) q[0];
rz(2.5671666) q[1];
sx q[1];
rz(-1.6822633) q[1];
sx q[1];
rz(-2.548545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7710073) q[0];
sx q[0];
rz(-0.86261504) q[0];
sx q[0];
rz(-0.38544356) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0306274) q[2];
sx q[2];
rz(-1.2048651) q[2];
sx q[2];
rz(2.1618787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30664602) q[1];
sx q[1];
rz(-3.1241841) q[1];
sx q[1];
rz(-2.2724292) q[1];
rz(-1.9022868) q[3];
sx q[3];
rz(-0.96486366) q[3];
sx q[3];
rz(-0.87406033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9888931) q[2];
sx q[2];
rz(-1.9579192) q[2];
sx q[2];
rz(-0.93842554) q[2];
rz(1.7484131) q[3];
sx q[3];
rz(-1.3104855) q[3];
sx q[3];
rz(-2.9032629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68168454) q[0];
sx q[0];
rz(-2.5354698) q[0];
sx q[0];
rz(1.348362) q[0];
rz(-2.0943191) q[1];
sx q[1];
rz(-0.92887639) q[1];
sx q[1];
rz(2.1705719) q[1];
rz(-2.1793096) q[2];
sx q[2];
rz(-3.0627927) q[2];
sx q[2];
rz(0.40712955) q[2];
rz(-1.3503475) q[3];
sx q[3];
rz(-1.7488283) q[3];
sx q[3];
rz(-1.8976952) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
