OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(-0.14578851) q[0];
sx q[0];
rz(1.2989651) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(4.4823449) q[1];
sx q[1];
rz(9.525099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31831384) q[0];
sx q[0];
rz(-0.66368503) q[0];
sx q[0];
rz(-0.084245988) q[0];
x q[1];
rz(1.9062564) q[2];
sx q[2];
rz(-2.0192207) q[2];
sx q[2];
rz(-1.8775196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70182204) q[1];
sx q[1];
rz(-1.3622829) q[1];
sx q[1];
rz(0.28065248) q[1];
rz(1.2572977) q[3];
sx q[3];
rz(-2.1380205) q[3];
sx q[3];
rz(2.6424266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7413062) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(-2.3347704) q[2];
rz(-0.72871366) q[3];
sx q[3];
rz(-2.3829134) q[3];
sx q[3];
rz(1.9369283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62419409) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(0.30701315) q[0];
rz(-1.2767876) q[1];
sx q[1];
rz(-1.1390319) q[1];
sx q[1];
rz(2.0397287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0372504) q[0];
sx q[0];
rz(-0.58339632) q[0];
sx q[0];
rz(-1.8895545) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6346143) q[2];
sx q[2];
rz(-1.7498657) q[2];
sx q[2];
rz(0.049953559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15436048) q[1];
sx q[1];
rz(-1.7061966) q[1];
sx q[1];
rz(1.7495278) q[1];
rz(-pi) q[2];
rz(0.73429196) q[3];
sx q[3];
rz(-2.2059388) q[3];
sx q[3];
rz(-1.7226302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1133984) q[2];
sx q[2];
rz(-0.58176175) q[2];
sx q[2];
rz(0.096435189) q[2];
rz(-0.15245572) q[3];
sx q[3];
rz(-1.6343445) q[3];
sx q[3];
rz(0.69795394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089712791) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(0.40019792) q[0];
rz(1.6290889) q[1];
sx q[1];
rz(-2.9632603) q[1];
sx q[1];
rz(2.8834744) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011964762) q[0];
sx q[0];
rz(-1.8045549) q[0];
sx q[0];
rz(-0.14272228) q[0];
rz(-1.6506972) q[2];
sx q[2];
rz(-1.4712442) q[2];
sx q[2];
rz(0.71074206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5096881) q[1];
sx q[1];
rz(-3.1277165) q[1];
sx q[1];
rz(-0.3664151) q[1];
rz(-pi) q[2];
rz(0.24127023) q[3];
sx q[3];
rz(-1.7381958) q[3];
sx q[3];
rz(-0.77360717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.021598024) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(3.1094587) q[2];
rz(2.8739127) q[3];
sx q[3];
rz(-1.5175502) q[3];
sx q[3];
rz(-1.5599498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5699919) q[0];
sx q[0];
rz(-1.6331693) q[0];
sx q[0];
rz(-0.084240325) q[0];
rz(-0.035942297) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(-2.8004004) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0628898) q[0];
sx q[0];
rz(-2.4475532) q[0];
sx q[0];
rz(2.7258123) q[0];
rz(1.8282537) q[2];
sx q[2];
rz(-0.99960589) q[2];
sx q[2];
rz(-1.5549623) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.045121) q[1];
sx q[1];
rz(-2.0092138) q[1];
sx q[1];
rz(-2.3271266) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59935948) q[3];
sx q[3];
rz(-1.9512358) q[3];
sx q[3];
rz(-2.3886556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7146032) q[2];
sx q[2];
rz(-2.0689071) q[2];
sx q[2];
rz(-1.6061456) q[2];
rz(2.8793907) q[3];
sx q[3];
rz(-1.6180792) q[3];
sx q[3];
rz(-0.95482993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.75761211) q[0];
sx q[0];
rz(-2.7181427) q[0];
sx q[0];
rz(2.0641932) q[0];
rz(2.7491838) q[1];
sx q[1];
rz(-0.078374021) q[1];
sx q[1];
rz(1.0307301) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1218613) q[0];
sx q[0];
rz(-1.9872338) q[0];
sx q[0];
rz(0.58953489) q[0];
rz(-1.0723128) q[2];
sx q[2];
rz(-1.7297812) q[2];
sx q[2];
rz(2.0778401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3001316) q[1];
sx q[1];
rz(-2.8634954) q[1];
sx q[1];
rz(0.93021955) q[1];
rz(1.0593653) q[3];
sx q[3];
rz(-1.2687195) q[3];
sx q[3];
rz(1.6405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3173759) q[2];
sx q[2];
rz(-2.4874096) q[2];
sx q[2];
rz(2.3023494) q[2];
rz(-0.9134891) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.6827253) q[0];
sx q[0];
rz(-2.904628) q[0];
sx q[0];
rz(-1.6656026) q[0];
rz(-2.7492375) q[1];
sx q[1];
rz(-2.0456435) q[1];
sx q[1];
rz(2.5501693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77350835) q[0];
sx q[0];
rz(-0.71660935) q[0];
sx q[0];
rz(-2.1518097) q[0];
rz(-pi) q[1];
rz(0.27523454) q[2];
sx q[2];
rz(-1.396759) q[2];
sx q[2];
rz(1.4402207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.082669584) q[1];
sx q[1];
rz(-1.6367404) q[1];
sx q[1];
rz(-1.0643427) q[1];
rz(0.87402451) q[3];
sx q[3];
rz(-1.2748147) q[3];
sx q[3];
rz(-0.10914739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8412987) q[2];
sx q[2];
rz(-0.66606194) q[2];
sx q[2];
rz(-2.5692614) q[2];
rz(0.22219292) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7543024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(2.733316) q[0];
rz(2.4093742) q[1];
sx q[1];
rz(-3.0156942) q[1];
sx q[1];
rz(-0.29762038) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10935303) q[0];
sx q[0];
rz(-1.223757) q[0];
sx q[0];
rz(0.45022398) q[0];
rz(-1.5333129) q[2];
sx q[2];
rz(-2.6465694) q[2];
sx q[2];
rz(-1.2353473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0875594) q[1];
sx q[1];
rz(-0.39759025) q[1];
sx q[1];
rz(2.1000186) q[1];
rz(-pi) q[2];
rz(1.7209531) q[3];
sx q[3];
rz(-2.2689226) q[3];
sx q[3];
rz(-1.0696749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1709661) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(-0.69620281) q[2];
rz(-2.2512186) q[3];
sx q[3];
rz(-1.1768769) q[3];
sx q[3];
rz(-1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9625229) q[0];
sx q[0];
rz(-0.028554976) q[0];
sx q[0];
rz(-0.21275511) q[0];
rz(2.6720324) q[1];
sx q[1];
rz(-0.9649562) q[1];
sx q[1];
rz(-0.75417095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6109723) q[0];
sx q[0];
rz(-1.7585131) q[0];
sx q[0];
rz(-1.2769481) q[0];
rz(-pi) q[1];
rz(1.7751664) q[2];
sx q[2];
rz(-1.8897383) q[2];
sx q[2];
rz(-1.5103923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2060125) q[1];
sx q[1];
rz(-1.4569062) q[1];
sx q[1];
rz(-2.4920032) q[1];
rz(-pi) q[2];
rz(0.92408224) q[3];
sx q[3];
rz(-1.0093186) q[3];
sx q[3];
rz(1.1738861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.135005) q[2];
sx q[2];
rz(-0.90234119) q[2];
sx q[2];
rz(-2.3593486) q[2];
rz(1.6953281) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(-2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049147216) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(-0.96442047) q[0];
rz(1.8655221) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.5020348) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0014711) q[0];
sx q[0];
rz(-1.7709416) q[0];
sx q[0];
rz(1.3339304) q[0];
rz(1.1310857) q[2];
sx q[2];
rz(-0.88958101) q[2];
sx q[2];
rz(2.0063673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0037476841) q[1];
sx q[1];
rz(-1.5478351) q[1];
sx q[1];
rz(-2.0883043) q[1];
rz(-pi) q[2];
rz(3.0739932) q[3];
sx q[3];
rz(-1.7177594) q[3];
sx q[3];
rz(2.802668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8883349) q[2];
sx q[2];
rz(-1.8586681) q[2];
sx q[2];
rz(0.9453195) q[2];
rz(0.82593289) q[3];
sx q[3];
rz(-1.632894) q[3];
sx q[3];
rz(-0.53502423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076040529) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(0.8031351) q[0];
rz(-1.5638634) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(2.8520083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591551) q[0];
sx q[0];
rz(-1.6763601) q[0];
sx q[0];
rz(-3.1292874) q[0];
x q[1];
rz(-2.824823) q[2];
sx q[2];
rz(-1.8837591) q[2];
sx q[2];
rz(0.53887689) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8378505) q[1];
sx q[1];
rz(-1.2134064) q[1];
sx q[1];
rz(-2.9010495) q[1];
x q[2];
rz(1.9042468) q[3];
sx q[3];
rz(-1.8419918) q[3];
sx q[3];
rz(-1.9681794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3760066) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(-0.96013367) q[2];
rz(-2.5586186) q[3];
sx q[3];
rz(-2.4834902) q[3];
sx q[3];
rz(-0.8647024) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4509907) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(-3.1008537) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(2.8340451) q[2];
sx q[2];
rz(-2.0219621) q[2];
sx q[2];
rz(-2.1297217) q[2];
rz(0.05091713) q[3];
sx q[3];
rz(-1.6418727) q[3];
sx q[3];
rz(1.4486936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
