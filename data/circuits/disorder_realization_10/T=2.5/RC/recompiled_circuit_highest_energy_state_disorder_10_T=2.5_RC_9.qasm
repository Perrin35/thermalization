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
rz(-0.69775692) q[0];
sx q[0];
rz(-1.193576) q[0];
sx q[0];
rz(-0.20456631) q[0];
rz(2.3808631) q[1];
sx q[1];
rz(-1.7525571) q[1];
sx q[1];
rz(-1.3995481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4142864) q[0];
sx q[0];
rz(-2.901863) q[0];
sx q[0];
rz(1.9637462) q[0];
rz(-2.2967161) q[2];
sx q[2];
rz(-2.5683218) q[2];
sx q[2];
rz(2.5115761) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1408994) q[1];
sx q[1];
rz(-1.4834036) q[1];
sx q[1];
rz(-2.7192975) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3035965) q[3];
sx q[3];
rz(-2.2624216) q[3];
sx q[3];
rz(1.8593781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7064887) q[2];
sx q[2];
rz(-3.1092643) q[2];
sx q[2];
rz(-2.6522563) q[2];
rz(2.1183744) q[3];
sx q[3];
rz(-0.018298572) q[3];
sx q[3];
rz(-2.0160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045687549) q[0];
sx q[0];
rz(-0.65653312) q[0];
sx q[0];
rz(-0.80279654) q[0];
rz(-0.071391694) q[1];
sx q[1];
rz(-2.8749021) q[1];
sx q[1];
rz(0.057770483) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5934385) q[0];
sx q[0];
rz(-1.8280941) q[0];
sx q[0];
rz(-2.1307039) q[0];
x q[1];
rz(2.7947756) q[2];
sx q[2];
rz(-0.3036193) q[2];
sx q[2];
rz(-0.70250073) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7747468) q[1];
sx q[1];
rz(-0.54301942) q[1];
sx q[1];
rz(1.1096891) q[1];
rz(-pi) q[2];
rz(2.9936419) q[3];
sx q[3];
rz(-0.86455621) q[3];
sx q[3];
rz(2.1220292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1959261) q[2];
sx q[2];
rz(-1.095093) q[2];
sx q[2];
rz(1.8460974) q[2];
rz(-0.96674353) q[3];
sx q[3];
rz(-0.77015489) q[3];
sx q[3];
rz(-0.75743341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.029723786) q[0];
sx q[0];
rz(-1.3628549) q[0];
sx q[0];
rz(-1.7080074) q[0];
rz(-3.0729821) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(-2.5624018) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5369956) q[0];
sx q[0];
rz(-1.662435) q[0];
sx q[0];
rz(-3.0834404) q[0];
x q[1];
rz(-1.5689108) q[2];
sx q[2];
rz(-2.0465188) q[2];
sx q[2];
rz(1.204551) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2229742) q[1];
sx q[1];
rz(-0.22906216) q[1];
sx q[1];
rz(0.15840662) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0924358) q[3];
sx q[3];
rz(-1.98171) q[3];
sx q[3];
rz(2.876407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.27089831) q[2];
sx q[2];
rz(-1.2535932) q[2];
sx q[2];
rz(2.9581621) q[2];
rz(-2.3373248) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(2.6075294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19876984) q[0];
sx q[0];
rz(-3.0267921) q[0];
sx q[0];
rz(0.59952366) q[0];
rz(-0.41246688) q[1];
sx q[1];
rz(-0.020655276) q[1];
sx q[1];
rz(-2.1309158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3314047) q[0];
sx q[0];
rz(-1.4443026) q[0];
sx q[0];
rz(2.7888915) q[0];
rz(-pi) q[1];
rz(-1.640219) q[2];
sx q[2];
rz(-0.63328082) q[2];
sx q[2];
rz(-0.4236003) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0914826) q[1];
sx q[1];
rz(-2.1569139) q[1];
sx q[1];
rz(0.32663235) q[1];
rz(-pi) q[2];
rz(-2.9948576) q[3];
sx q[3];
rz(-0.76125604) q[3];
sx q[3];
rz(3.0748607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3707054) q[2];
sx q[2];
rz(-0.34660307) q[2];
sx q[2];
rz(2.8240805) q[2];
rz(2.5192449) q[3];
sx q[3];
rz(-0.93572664) q[3];
sx q[3];
rz(0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7790826) q[0];
sx q[0];
rz(-2.2267987) q[0];
sx q[0];
rz(-2.1146178) q[0];
rz(0.55038553) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(-1.9245573) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7893697) q[0];
sx q[0];
rz(-2.2055211) q[0];
sx q[0];
rz(2.288649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53834277) q[2];
sx q[2];
rz(-1.6377047) q[2];
sx q[2];
rz(2.4198893) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74351013) q[1];
sx q[1];
rz(-2.4226502) q[1];
sx q[1];
rz(1.8330857) q[1];
rz(-2.2812165) q[3];
sx q[3];
rz(-1.2902707) q[3];
sx q[3];
rz(-0.70395148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43651849) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(2.3684033) q[2];
rz(-3.0298722) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(2.2022061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47950995) q[0];
sx q[0];
rz(-2.918512) q[0];
sx q[0];
rz(2.7146085) q[0];
rz(-0.93049479) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(0.46447909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31373608) q[0];
sx q[0];
rz(-1.858874) q[0];
sx q[0];
rz(-1.366607) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1709178) q[2];
sx q[2];
rz(-1.8541186) q[2];
sx q[2];
rz(1.4221869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.528094) q[1];
sx q[1];
rz(-1.6069176) q[1];
sx q[1];
rz(0.36904676) q[1];
rz(-pi) q[2];
rz(1.6397912) q[3];
sx q[3];
rz(-2.7164) q[3];
sx q[3];
rz(1.0591648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6102607) q[2];
sx q[2];
rz(-1.8208296) q[2];
sx q[2];
rz(-0.28826928) q[2];
rz(1.0432976) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(2.402795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6601324) q[0];
sx q[0];
rz(-1.8539424) q[0];
sx q[0];
rz(-2.3840391) q[0];
rz(0.07269147) q[1];
sx q[1];
rz(-0.025765954) q[1];
sx q[1];
rz(3.0933948) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.165503) q[0];
sx q[0];
rz(-0.58506008) q[0];
sx q[0];
rz(-1.6242322) q[0];
rz(-pi) q[1];
rz(-2.169007) q[2];
sx q[2];
rz(-1.3763104) q[2];
sx q[2];
rz(-2.011428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1485702) q[1];
sx q[1];
rz(-0.76495586) q[1];
sx q[1];
rz(-2.2123442) q[1];
rz(-pi) q[2];
rz(-0.052560135) q[3];
sx q[3];
rz(-2.0400999) q[3];
sx q[3];
rz(-1.8183501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4797719) q[2];
sx q[2];
rz(-1.6419819) q[2];
sx q[2];
rz(-3.0721967) q[2];
rz(1.5540468) q[3];
sx q[3];
rz(-0.78726751) q[3];
sx q[3];
rz(2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.7540392) q[0];
sx q[0];
rz(-1.0559005) q[0];
sx q[0];
rz(-1.7283424) q[0];
rz(-2.3130401) q[1];
sx q[1];
rz(-0.041752432) q[1];
sx q[1];
rz(-2.5989596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56930243) q[0];
sx q[0];
rz(-2.9539032) q[0];
sx q[0];
rz(-2.45298) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9285525) q[2];
sx q[2];
rz(-1.2160436) q[2];
sx q[2];
rz(-0.32217978) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9787776) q[1];
sx q[1];
rz(-0.12773027) q[1];
sx q[1];
rz(-2.7055897) q[1];
x q[2];
rz(-1.1922791) q[3];
sx q[3];
rz(-1.1633486) q[3];
sx q[3];
rz(0.70127869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98216206) q[2];
sx q[2];
rz(-0.40951481) q[2];
sx q[2];
rz(2.8816667) q[2];
rz(-2.121117) q[3];
sx q[3];
rz(-2.8838938) q[3];
sx q[3];
rz(-2.3475588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644153) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(1.6625241) q[0];
rz(1.5326477) q[1];
sx q[1];
rz(-1.0391935) q[1];
sx q[1];
rz(0.74302465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6830269) q[0];
sx q[0];
rz(-1.6855778) q[0];
sx q[0];
rz(1.0558788) q[0];
x q[1];
rz(1.0270306) q[2];
sx q[2];
rz(-1.2390422) q[2];
sx q[2];
rz(2.5912655) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1260657) q[1];
sx q[1];
rz(-3.0558476) q[1];
sx q[1];
rz(-2.2669906) q[1];
x q[2];
rz(3.0366692) q[3];
sx q[3];
rz(-1.6015953) q[3];
sx q[3];
rz(-1.2237751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6741901) q[2];
sx q[2];
rz(-0.82618606) q[2];
sx q[2];
rz(-0.80545938) q[2];
rz(-1.6921267) q[3];
sx q[3];
rz(-1.9129246) q[3];
sx q[3];
rz(-0.8031556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8886803) q[0];
sx q[0];
rz(-0.54179931) q[0];
sx q[0];
rz(-2.3895277) q[0];
rz(-1.1532785) q[1];
sx q[1];
rz(-2.2572932) q[1];
sx q[1];
rz(-2.8582252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3842938) q[0];
sx q[0];
rz(-1.481825) q[0];
sx q[0];
rz(-2.6020537) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4128017) q[2];
sx q[2];
rz(-0.758095) q[2];
sx q[2];
rz(-1.5811063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9787406) q[1];
sx q[1];
rz(-1.2417485) q[1];
sx q[1];
rz(0.63905893) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9612736) q[3];
sx q[3];
rz(-2.1877398) q[3];
sx q[3];
rz(2.4882567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.89455426) q[2];
sx q[2];
rz(-3.0587695) q[2];
sx q[2];
rz(-1.7029597) q[2];
rz(0.29397193) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(-2.1132052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033584874) q[0];
sx q[0];
rz(-1.7243732) q[0];
sx q[0];
rz(1.617817) q[0];
rz(-0.53957466) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(1.9488967) q[2];
sx q[2];
rz(-0.38496502) q[2];
sx q[2];
rz(-2.5733583) q[2];
rz(-1.2992819) q[3];
sx q[3];
rz(-1.9403024) q[3];
sx q[3];
rz(-1.3997146) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
