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
rz(-2.3186853) q[0];
sx q[0];
rz(-2.7828001) q[0];
sx q[0];
rz(-0.86831492) q[0];
rz(1.7262285) q[1];
sx q[1];
rz(-2.0077029) q[1];
sx q[1];
rz(0.99376065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9159587) q[0];
sx q[0];
rz(-1.6835815) q[0];
sx q[0];
rz(-2.7382572) q[0];
x q[1];
rz(-2.6401071) q[2];
sx q[2];
rz(-1.9619313) q[2];
sx q[2];
rz(2.8361965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3512416) q[1];
sx q[1];
rz(-2.0000441) q[1];
sx q[1];
rz(2.0584978) q[1];
x q[2];
rz(0.24561974) q[3];
sx q[3];
rz(-1.1455451) q[3];
sx q[3];
rz(-0.58477816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.913784) q[2];
sx q[2];
rz(-1.8081534) q[2];
sx q[2];
rz(-0.58161962) q[2];
rz(2.4070814) q[3];
sx q[3];
rz(-1.4903277) q[3];
sx q[3];
rz(2.7233126) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89473474) q[0];
sx q[0];
rz(-1.3792091) q[0];
sx q[0];
rz(2.6439164) q[0];
rz(1.0408164) q[1];
sx q[1];
rz(-2.7383995) q[1];
sx q[1];
rz(-1.9042447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055641551) q[0];
sx q[0];
rz(-2.4179672) q[0];
sx q[0];
rz(-0.063972278) q[0];
rz(1.2096268) q[2];
sx q[2];
rz(-1.7470844) q[2];
sx q[2];
rz(0.90106264) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.1488545) q[1];
sx q[1];
rz(-1.142456) q[1];
sx q[1];
rz(-0.15495877) q[1];
x q[2];
rz(1.353831) q[3];
sx q[3];
rz(-1.8304018) q[3];
sx q[3];
rz(-1.6093773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70025468) q[2];
sx q[2];
rz(-0.73107084) q[2];
sx q[2];
rz(2.8311484) q[2];
rz(-1.1566409) q[3];
sx q[3];
rz(-1.1247331) q[3];
sx q[3];
rz(-1.9748851) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0079086) q[0];
sx q[0];
rz(-2.5569361) q[0];
sx q[0];
rz(0.34580082) q[0];
rz(2.8969104) q[1];
sx q[1];
rz(-2.209765) q[1];
sx q[1];
rz(-1.4366879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61621746) q[0];
sx q[0];
rz(-1.8843009) q[0];
sx q[0];
rz(2.0140735) q[0];
rz(1.9635003) q[2];
sx q[2];
rz(-1.0356734) q[2];
sx q[2];
rz(1.2924043) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2785221) q[1];
sx q[1];
rz(-0.82048847) q[1];
sx q[1];
rz(0.29149518) q[1];
x q[2];
rz(2.6639943) q[3];
sx q[3];
rz(-1.2354038) q[3];
sx q[3];
rz(-0.23303495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1338542) q[2];
sx q[2];
rz(-1.1275007) q[2];
sx q[2];
rz(-0.63068843) q[2];
rz(-2.808029) q[3];
sx q[3];
rz(-2.1400698) q[3];
sx q[3];
rz(1.001531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075060464) q[0];
sx q[0];
rz(-2.1933031) q[0];
sx q[0];
rz(1.7370268) q[0];
rz(0.36918494) q[1];
sx q[1];
rz(-1.4063947) q[1];
sx q[1];
rz(-0.085748347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.551683) q[0];
sx q[0];
rz(-2.6250408) q[0];
sx q[0];
rz(0.63802436) q[0];
x q[1];
rz(0.71620415) q[2];
sx q[2];
rz(-2.5599481) q[2];
sx q[2];
rz(-3.1379791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4300116) q[1];
sx q[1];
rz(-2.0329355) q[1];
sx q[1];
rz(0.26600809) q[1];
x q[2];
rz(-0.39601456) q[3];
sx q[3];
rz(-1.2761444) q[3];
sx q[3];
rz(-2.9673607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3021476) q[2];
sx q[2];
rz(-1.2371233) q[2];
sx q[2];
rz(-2.6117924) q[2];
rz(-3.1032622) q[3];
sx q[3];
rz(-2.4119792) q[3];
sx q[3];
rz(-1.5420325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9180651) q[0];
sx q[0];
rz(-0.46850884) q[0];
sx q[0];
rz(1.202762) q[0];
rz(-0.27944061) q[1];
sx q[1];
rz(-1.025082) q[1];
sx q[1];
rz(0.7739982) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5764172) q[0];
sx q[0];
rz(-0.87067184) q[0];
sx q[0];
rz(-2.7267674) q[0];
rz(-pi) q[1];
rz(-1.4882795) q[2];
sx q[2];
rz(-1.8029717) q[2];
sx q[2];
rz(1.9407995) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8355636) q[1];
sx q[1];
rz(-1.4317498) q[1];
sx q[1];
rz(-3.1285888) q[1];
rz(-pi) q[2];
rz(-1.5660902) q[3];
sx q[3];
rz(-2.60618) q[3];
sx q[3];
rz(-2.0369867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0224132) q[2];
sx q[2];
rz(-3.000562) q[2];
sx q[2];
rz(-2.5775583) q[2];
rz(2.4885079) q[3];
sx q[3];
rz(-1.1354732) q[3];
sx q[3];
rz(-0.45026067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7804724) q[0];
sx q[0];
rz(-0.55279624) q[0];
sx q[0];
rz(0.64315382) q[0];
rz(-1.1910575) q[1];
sx q[1];
rz(-1.4570313) q[1];
sx q[1];
rz(-2.4868884) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1977256) q[0];
sx q[0];
rz(-1.7392613) q[0];
sx q[0];
rz(-1.6351624) q[0];
x q[1];
rz(-3.0949101) q[2];
sx q[2];
rz(-0.68447036) q[2];
sx q[2];
rz(-0.6990664) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1015849) q[1];
sx q[1];
rz(-1.5052649) q[1];
sx q[1];
rz(1.0270018) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1804462) q[3];
sx q[3];
rz(-1.6359513) q[3];
sx q[3];
rz(2.3303243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5074629) q[2];
sx q[2];
rz(-2.7644988) q[2];
sx q[2];
rz(1.6667574) q[2];
rz(2.6569488) q[3];
sx q[3];
rz(-2.1438997) q[3];
sx q[3];
rz(-1.8528329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.5019048) q[0];
sx q[0];
rz(-0.26308331) q[0];
sx q[0];
rz(0.68341533) q[0];
rz(-3.0112093) q[1];
sx q[1];
rz(-1.5510635) q[1];
sx q[1];
rz(0.032616671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17856971) q[0];
sx q[0];
rz(-1.7462329) q[0];
sx q[0];
rz(-1.001872) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24222272) q[2];
sx q[2];
rz(-0.70021473) q[2];
sx q[2];
rz(-1.2716573) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4144145) q[1];
sx q[1];
rz(-1.7373996) q[1];
sx q[1];
rz(1.3018621) q[1];
rz(-1.2682876) q[3];
sx q[3];
rz(-2.4166346) q[3];
sx q[3];
rz(2.8444949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7183097) q[2];
sx q[2];
rz(-2.0330567) q[2];
sx q[2];
rz(0.56524593) q[2];
rz(1.7806753) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(-0.81418532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77124202) q[0];
sx q[0];
rz(-3.0841565) q[0];
sx q[0];
rz(-2.8420319) q[0];
rz(1.4467422) q[1];
sx q[1];
rz(-2.2694777) q[1];
sx q[1];
rz(2.1902693) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0977391) q[0];
sx q[0];
rz(-1.7758177) q[0];
sx q[0];
rz(-0.033173843) q[0];
x q[1];
rz(3.1323593) q[2];
sx q[2];
rz(-1.7823185) q[2];
sx q[2];
rz(-0.68995014) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48288667) q[1];
sx q[1];
rz(-1.1466007) q[1];
sx q[1];
rz(-1.5931604) q[1];
x q[2];
rz(-1.8325915) q[3];
sx q[3];
rz(-1.6556532) q[3];
sx q[3];
rz(3.121162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1424554) q[2];
sx q[2];
rz(-0.74644011) q[2];
sx q[2];
rz(2.8738521) q[2];
rz(1.0033876) q[3];
sx q[3];
rz(-0.95971862) q[3];
sx q[3];
rz(-0.80823922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830911) q[0];
sx q[0];
rz(-0.49810228) q[0];
sx q[0];
rz(-2.1844693) q[0];
rz(-0.3793017) q[1];
sx q[1];
rz(-2.7298268) q[1];
sx q[1];
rz(-0.1651102) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50277621) q[0];
sx q[0];
rz(-0.90327016) q[0];
sx q[0];
rz(1.1291885) q[0];
x q[1];
rz(1.0983989) q[2];
sx q[2];
rz(-1.6159) q[2];
sx q[2];
rz(-1.2055604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2735426) q[1];
sx q[1];
rz(-0.99813491) q[1];
sx q[1];
rz(-0.21577253) q[1];
x q[2];
rz(0.96747193) q[3];
sx q[3];
rz(-1.3730197) q[3];
sx q[3];
rz(-2.8654049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86753201) q[2];
sx q[2];
rz(-1.2313077) q[2];
sx q[2];
rz(0.77862281) q[2];
rz(2.8042931) q[3];
sx q[3];
rz(-2.1179347) q[3];
sx q[3];
rz(-1.5571099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839277) q[0];
sx q[0];
rz(-2.0511257) q[0];
sx q[0];
rz(2.0528059) q[0];
rz(-2.7133443) q[1];
sx q[1];
rz(-1.0232404) q[1];
sx q[1];
rz(-0.64839378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.245609) q[0];
sx q[0];
rz(-1.147384) q[0];
sx q[0];
rz(-1.8603252) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9087547) q[2];
sx q[2];
rz(-1.612886) q[2];
sx q[2];
rz(1.011285) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32247816) q[1];
sx q[1];
rz(-1.6129061) q[1];
sx q[1];
rz(-2.4933715) q[1];
rz(-pi) q[2];
rz(-1.4565868) q[3];
sx q[3];
rz(-2.1314959) q[3];
sx q[3];
rz(1.4637092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8757561) q[2];
sx q[2];
rz(-2.0529604) q[2];
sx q[2];
rz(2.9618373) q[2];
rz(1.1579375) q[3];
sx q[3];
rz(-1.4561184) q[3];
sx q[3];
rz(1.2402844) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801878) q[0];
sx q[0];
rz(-2.9869933) q[0];
sx q[0];
rz(-2.3416478) q[0];
rz(1.9752621) q[1];
sx q[1];
rz(-0.98465289) q[1];
sx q[1];
rz(-2.226895) q[1];
rz(-2.2565319) q[2];
sx q[2];
rz(-0.37715465) q[2];
sx q[2];
rz(-2.5436795) q[2];
rz(2.8020482) q[3];
sx q[3];
rz(-1.3673269) q[3];
sx q[3];
rz(0.23585503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
