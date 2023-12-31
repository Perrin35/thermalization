OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(-1.3347081) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53544331) q[0];
sx q[0];
rz(-0.9978928) q[0];
sx q[0];
rz(1.2826305) q[0];
rz(-1.0078148) q[2];
sx q[2];
rz(-1.286176) q[2];
sx q[2];
rz(2.9993338) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18260278) q[1];
sx q[1];
rz(-0.62287736) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(-0.39335143) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(0.8918744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(-0.33256724) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(-1.8012841) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48288229) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20992499) q[0];
sx q[0];
rz(-1.4405662) q[0];
sx q[0];
rz(1.7658712) q[0];
rz(-pi) q[1];
rz(-0.52829929) q[2];
sx q[2];
rz(-1.0353147) q[2];
sx q[2];
rz(-2.2250125) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43286846) q[1];
sx q[1];
rz(-1.8808865) q[1];
sx q[1];
rz(1.6354145) q[1];
rz(-pi) q[2];
rz(-0.56224058) q[3];
sx q[3];
rz(-1.3573109) q[3];
sx q[3];
rz(-1.436304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(2.7745461) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(0.095741622) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.384882) q[0];
sx q[0];
rz(-1.7236992) q[0];
sx q[0];
rz(1.3923313) q[0];
x q[1];
rz(0.50600608) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(2.4207052) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.15624554) q[1];
sx q[1];
rz(-2.3932082) q[1];
sx q[1];
rz(2.2845539) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0843094) q[3];
sx q[3];
rz(-0.55192845) q[3];
sx q[3];
rz(-2.3416167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(2.3068008) q[2];
rz(0.21162027) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(-0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-2.7950177) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(-1.0338763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047065145) q[0];
sx q[0];
rz(-2.3405502) q[0];
sx q[0];
rz(-0.69181504) q[0];
rz(-0.526555) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(-1.9499792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.911072) q[1];
sx q[1];
rz(-1.5057179) q[1];
sx q[1];
rz(-1.7622403) q[1];
x q[2];
rz(2.2771308) q[3];
sx q[3];
rz(-2.7221788) q[3];
sx q[3];
rz(-2.1692587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(-2.7205617) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(1.746009) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(-2.5674852) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610883) q[0];
sx q[0];
rz(-1.442712) q[0];
sx q[0];
rz(-2.1696521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5482043) q[2];
sx q[2];
rz(-2.2859757) q[2];
sx q[2];
rz(1.9695645) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9628323) q[1];
sx q[1];
rz(-1.8492336) q[1];
sx q[1];
rz(0.36414418) q[1];
rz(-1.3247213) q[3];
sx q[3];
rz(-2.4613791) q[3];
sx q[3];
rz(-2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.5138907) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(-0.73227698) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9376611) q[0];
sx q[0];
rz(-1.7294149) q[0];
sx q[0];
rz(-2.926814) q[0];
rz(-2.4562624) q[2];
sx q[2];
rz(-2.2517859) q[2];
sx q[2];
rz(-2.6076917) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4477168) q[1];
sx q[1];
rz(-2.4442721) q[1];
sx q[1];
rz(-0.35481528) q[1];
rz(-pi) q[2];
rz(-1.3527649) q[3];
sx q[3];
rz(-1.0112959) q[3];
sx q[3];
rz(2.7263209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(-1.1614655) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(2.9851595) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(0.77004534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55243385) q[0];
sx q[0];
rz(-1.7434412) q[0];
sx q[0];
rz(2.5347559) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1595575) q[2];
sx q[2];
rz(-1.9141478) q[2];
sx q[2];
rz(2.8729168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8416482) q[1];
sx q[1];
rz(-2.3023459) q[1];
sx q[1];
rz(-2.4323835) q[1];
x q[2];
rz(1.2392063) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5148233) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-0.36639211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241656) q[0];
sx q[0];
rz(-0.96091849) q[0];
sx q[0];
rz(2.2141371) q[0];
rz(1.1364469) q[2];
sx q[2];
rz(-2.1546116) q[2];
sx q[2];
rz(0.39603147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3185127) q[1];
sx q[1];
rz(-1.4050583) q[1];
sx q[1];
rz(-0.70648273) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1226235) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(-3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(-2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-0.88395399) q[0];
rz(-2.4813095) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(-2.3502137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26830772) q[0];
sx q[0];
rz(-0.76457667) q[0];
sx q[0];
rz(1.1029878) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0558526) q[2];
sx q[2];
rz(-1.9888478) q[2];
sx q[2];
rz(-0.56318356) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1014175) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(-2.0324213) q[1];
rz(-2.1095554) q[3];
sx q[3];
rz(-1.4269392) q[3];
sx q[3];
rz(2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(2.1203314) q[2];
rz(2.7630473) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(2.9027477) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(1.7609319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35683435) q[0];
sx q[0];
rz(-2.1958302) q[0];
sx q[0];
rz(-2.9805095) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0062749) q[2];
sx q[2];
rz(-1.5473817) q[2];
sx q[2];
rz(-1.5362816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7271991) q[1];
sx q[1];
rz(-2.6135923) q[1];
sx q[1];
rz(-0.90331932) q[1];
rz(-pi) q[2];
rz(-1.4019743) q[3];
sx q[3];
rz(-1.3843378) q[3];
sx q[3];
rz(-2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-0.33622462) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1375785) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(0.54429383) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(2.8019194) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(1.929677) q[3];
sx q[3];
rz(-2.5419895) q[3];
sx q[3];
rz(0.06751577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
