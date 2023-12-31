OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(-2.3242216) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(-1.2004381) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6450206) q[0];
sx q[0];
rz(-1.3294157) q[0];
sx q[0];
rz(0.42005959) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1098233) q[2];
sx q[2];
rz(-1.6845778) q[2];
sx q[2];
rz(-0.10345085) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8564954) q[1];
sx q[1];
rz(-2.3772117) q[1];
sx q[1];
rz(-2.7544423) q[1];
x q[2];
rz(0.31906268) q[3];
sx q[3];
rz(-2.3205119) q[3];
sx q[3];
rz(1.0533489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1203221) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(1.583064) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(-0.054071991) q[0];
rz(-1.1955098) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(-0.53584677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065952259) q[0];
sx q[0];
rz(-1.6951121) q[0];
sx q[0];
rz(-2.1527704) q[0];
rz(-1.8230121) q[2];
sx q[2];
rz(-2.2644342) q[2];
sx q[2];
rz(1.1064305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51551137) q[1];
sx q[1];
rz(-1.9768081) q[1];
sx q[1];
rz(-1.2938234) q[1];
rz(-3.1090118) q[3];
sx q[3];
rz(-1.7136646) q[3];
sx q[3];
rz(-0.25861614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(2.1832441) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(-0.025807468) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8255071) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(-3.120963) q[0];
rz(1.843156) q[2];
sx q[2];
rz(-0.32264999) q[2];
sx q[2];
rz(1.1622365) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64297134) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(0.69795124) q[1];
rz(-pi) q[2];
rz(-3.0098626) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(-0.47117885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(2.5562111) q[2];
rz(2.9600926) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(1.5649149) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(-2.9649819) q[0];
rz(0.88090849) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(0.53612971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49318424) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(1.0113082) q[0];
rz(-pi) q[1];
rz(-1.5136112) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(-2.016071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94169468) q[1];
sx q[1];
rz(-0.75921339) q[1];
sx q[1];
rz(-0.56337507) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36066182) q[3];
sx q[3];
rz(-2.4324527) q[3];
sx q[3];
rz(-0.053645596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(2.0969351) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.9064643) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(1.0513603) q[0];
rz(-1.4936739) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(-0.043118127) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2152104) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(-0.2290639) q[0];
rz(0.28508614) q[2];
sx q[2];
rz(-1.0694155) q[2];
sx q[2];
rz(-1.3146871) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9777898) q[1];
sx q[1];
rz(-1.3063352) q[1];
sx q[1];
rz(2.8603641) q[1];
rz(-2.99302) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(-3.1032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0126426) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(-0.35219231) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5181638) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(-2.39134) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(1.0587143) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48460618) q[0];
sx q[0];
rz(-1.3625506) q[0];
sx q[0];
rz(-2.7013742) q[0];
rz(-pi) q[1];
rz(0.58745678) q[2];
sx q[2];
rz(-2.537478) q[2];
sx q[2];
rz(0.47746745) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42029542) q[1];
sx q[1];
rz(-1.2696206) q[1];
sx q[1];
rz(-2.9073614) q[1];
x q[2];
rz(1.3771463) q[3];
sx q[3];
rz(-2.2691257) q[3];
sx q[3];
rz(2.8805672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(-1.2188101) q[2];
rz(1.9865215) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999775) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(2.2999433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8182897) q[0];
sx q[0];
rz(-1.0536195) q[0];
sx q[0];
rz(1.3099758) q[0];
x q[1];
rz(-1.3782578) q[2];
sx q[2];
rz(-2.1888869) q[2];
sx q[2];
rz(2.3877909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0824273) q[1];
sx q[1];
rz(-1.4945684) q[1];
sx q[1];
rz(-0.04870292) q[1];
x q[2];
rz(0.82960415) q[3];
sx q[3];
rz(-1.4083574) q[3];
sx q[3];
rz(3.0138569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.001361751) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(-0.049953071) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426303) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(-1.4021953) q[0];
rz(0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1091052) q[0];
sx q[0];
rz(-1.0787449) q[0];
sx q[0];
rz(1.0743272) q[0];
x q[1];
rz(-1.278644) q[2];
sx q[2];
rz(-2.7610965) q[2];
sx q[2];
rz(-0.7064864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8199181) q[1];
sx q[1];
rz(-0.770861) q[1];
sx q[1];
rz(-0.8701156) q[1];
x q[2];
rz(2.8475259) q[3];
sx q[3];
rz(-1.9293647) q[3];
sx q[3];
rz(-0.29953526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(-2.4677094) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8326571) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(2.5667403) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(-0.94690698) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6553584) q[0];
sx q[0];
rz(-1.8706733) q[0];
sx q[0];
rz(-2.8006058) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5309179) q[2];
sx q[2];
rz(-1.6677742) q[2];
sx q[2];
rz(3.1181042) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0526035) q[1];
sx q[1];
rz(-0.34480428) q[1];
sx q[1];
rz(0.25119541) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1831207) q[3];
sx q[3];
rz(-2.6689853) q[3];
sx q[3];
rz(1.8567059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(-0.71436626) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0666075) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(-2.4865436) q[0];
rz(0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(-2.4972829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9172168) q[0];
sx q[0];
rz(-1.3972358) q[0];
sx q[0];
rz(-1.6266842) q[0];
x q[1];
rz(0.23004736) q[2];
sx q[2];
rz(-0.86798475) q[2];
sx q[2];
rz(2.7237797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88313738) q[1];
sx q[1];
rz(-1.2346134) q[1];
sx q[1];
rz(-1.2698445) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82018567) q[3];
sx q[3];
rz(-2.5327442) q[3];
sx q[3];
rz(-1.604515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(-1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(-3.0145666) q[2];
sx q[2];
rz(-0.95022485) q[2];
sx q[2];
rz(-2.7765204) q[2];
rz(2.7502144) q[3];
sx q[3];
rz(-2.2073675) q[3];
sx q[3];
rz(-1.7195306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
