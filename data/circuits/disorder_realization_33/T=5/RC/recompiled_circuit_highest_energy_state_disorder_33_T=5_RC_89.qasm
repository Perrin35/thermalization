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
rz(1.1627816) q[0];
sx q[0];
rz(-0.94036189) q[0];
sx q[0];
rz(-2.9475687) q[0];
rz(2.5972875) q[1];
sx q[1];
rz(-1.5736009) q[1];
sx q[1];
rz(0.1967217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5797884) q[0];
sx q[0];
rz(-1.6937979) q[0];
sx q[0];
rz(-2.1462026) q[0];
x q[1];
rz(2.9241184) q[2];
sx q[2];
rz(-1.5921235) q[2];
sx q[2];
rz(0.77922098) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87952033) q[1];
sx q[1];
rz(-1.2517057) q[1];
sx q[1];
rz(-1.2882922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3338577) q[3];
sx q[3];
rz(-0.58962017) q[3];
sx q[3];
rz(-0.20649621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0611614) q[2];
sx q[2];
rz(-1.4208527) q[2];
sx q[2];
rz(-0.013193456) q[2];
rz(2.0773928) q[3];
sx q[3];
rz(-2.1049757) q[3];
sx q[3];
rz(-0.17087759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5203633) q[0];
sx q[0];
rz(-0.55084387) q[0];
sx q[0];
rz(-2.6589822) q[0];
rz(-2.6699325) q[1];
sx q[1];
rz(-0.99716798) q[1];
sx q[1];
rz(-2.4536123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9887138) q[0];
sx q[0];
rz(-1.7611302) q[0];
sx q[0];
rz(2.8318263) q[0];
rz(-1.2763763) q[2];
sx q[2];
rz(-2.3465119) q[2];
sx q[2];
rz(-2.9816422) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6370442) q[1];
sx q[1];
rz(-2.6682862) q[1];
sx q[1];
rz(3.0454703) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1766731) q[3];
sx q[3];
rz(-2.9218194) q[3];
sx q[3];
rz(-3.0434853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1241577) q[2];
sx q[2];
rz(-1.0729125) q[2];
sx q[2];
rz(-2.9166481) q[2];
rz(2.0624835) q[3];
sx q[3];
rz(-0.93828097) q[3];
sx q[3];
rz(-0.45043954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17063046) q[0];
sx q[0];
rz(-0.58929515) q[0];
sx q[0];
rz(1.851409) q[0];
rz(1.2782485) q[1];
sx q[1];
rz(-1.1561013) q[1];
sx q[1];
rz(-1.2826756) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0541924) q[0];
sx q[0];
rz(-2.4303959) q[0];
sx q[0];
rz(-1.0048701) q[0];
x q[1];
rz(-3.01802) q[2];
sx q[2];
rz(-1.4915492) q[2];
sx q[2];
rz(2.045897) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45047347) q[1];
sx q[1];
rz(-1.5509924) q[1];
sx q[1];
rz(3.1234689) q[1];
rz(-0.76126666) q[3];
sx q[3];
rz(-2.1652615) q[3];
sx q[3];
rz(0.52235421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67905417) q[2];
sx q[2];
rz(-1.6109698) q[2];
sx q[2];
rz(-0.47492096) q[2];
rz(-1.9654407) q[3];
sx q[3];
rz(-2.2852496) q[3];
sx q[3];
rz(-0.28158751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7684286) q[0];
sx q[0];
rz(-1.3685065) q[0];
sx q[0];
rz(0.1678634) q[0];
rz(-2.8236875) q[1];
sx q[1];
rz(-1.8524421) q[1];
sx q[1];
rz(2.1071404) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7307593) q[0];
sx q[0];
rz(-0.26450464) q[0];
sx q[0];
rz(1.4368617) q[0];
rz(-pi) q[1];
rz(1.9475161) q[2];
sx q[2];
rz(-0.92096201) q[2];
sx q[2];
rz(-1.5638994) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.381731) q[1];
sx q[1];
rz(-1.1118366) q[1];
sx q[1];
rz(-0.81736395) q[1];
x q[2];
rz(0.14334821) q[3];
sx q[3];
rz(-0.61666538) q[3];
sx q[3];
rz(0.35307742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0340524) q[2];
sx q[2];
rz(-0.62959051) q[2];
sx q[2];
rz(-0.26629392) q[2];
rz(-2.9722424) q[3];
sx q[3];
rz(-1.063187) q[3];
sx q[3];
rz(0.17899409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-2.2696445) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(-1.2931152) q[0];
rz(-3.0710908) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(2.3354796) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370395) q[0];
sx q[0];
rz(-1.4463639) q[0];
sx q[0];
rz(1.7794153) q[0];
x q[1];
rz(0.057439645) q[2];
sx q[2];
rz(-1.7195133) q[2];
sx q[2];
rz(1.1401389) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5841585) q[1];
sx q[1];
rz(-0.63156742) q[1];
sx q[1];
rz(2.8110792) q[1];
rz(-pi) q[2];
rz(1.7905037) q[3];
sx q[3];
rz(-1.0423805) q[3];
sx q[3];
rz(-2.3912969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9301794) q[2];
sx q[2];
rz(-1.4750007) q[2];
sx q[2];
rz(-2.831947) q[2];
rz(0.51112255) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(-2.3206594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9822134) q[0];
sx q[0];
rz(-0.62726227) q[0];
sx q[0];
rz(-0.43701592) q[0];
rz(0.45477319) q[1];
sx q[1];
rz(-0.79650703) q[1];
sx q[1];
rz(2.2589267) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012271492) q[0];
sx q[0];
rz(-2.8573749) q[0];
sx q[0];
rz(1.8924665) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1922311) q[2];
sx q[2];
rz(-1.6003216) q[2];
sx q[2];
rz(-1.6395417) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5458293) q[1];
sx q[1];
rz(-0.59844164) q[1];
sx q[1];
rz(-0.95495895) q[1];
x q[2];
rz(2.1737133) q[3];
sx q[3];
rz(-1.2025739) q[3];
sx q[3];
rz(-2.8944824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5060045) q[2];
sx q[2];
rz(-1.4281861) q[2];
sx q[2];
rz(-0.44090718) q[2];
rz(-2.0697557) q[3];
sx q[3];
rz(-0.25522885) q[3];
sx q[3];
rz(2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275682) q[0];
sx q[0];
rz(-1.9381645) q[0];
sx q[0];
rz(-1.8593651) q[0];
rz(-2.0558689) q[1];
sx q[1];
rz(-1.6174569) q[1];
sx q[1];
rz(-1.5132743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0659637) q[0];
sx q[0];
rz(-0.83978486) q[0];
sx q[0];
rz(-2.9828067) q[0];
rz(-pi) q[1];
rz(2.7658913) q[2];
sx q[2];
rz(-0.73601572) q[2];
sx q[2];
rz(-1.0003566) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5243064) q[1];
sx q[1];
rz(-0.77037659) q[1];
sx q[1];
rz(-1.0138033) q[1];
rz(0.75700883) q[3];
sx q[3];
rz(-1.8534213) q[3];
sx q[3];
rz(-2.4429537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6699803) q[2];
sx q[2];
rz(-0.70576224) q[2];
sx q[2];
rz(2.3213049) q[2];
rz(-0.39508501) q[3];
sx q[3];
rz(-2.3494651) q[3];
sx q[3];
rz(0.61409605) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87913269) q[0];
sx q[0];
rz(-1.6430055) q[0];
sx q[0];
rz(2.5347128) q[0];
rz(-0.34582368) q[1];
sx q[1];
rz(-1.9551829) q[1];
sx q[1];
rz(-0.80723673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5257634) q[0];
sx q[0];
rz(-0.35789546) q[0];
sx q[0];
rz(-1.8496022) q[0];
x q[1];
rz(-1.973602) q[2];
sx q[2];
rz(-1.323408) q[2];
sx q[2];
rz(-2.0217264) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9620558) q[1];
sx q[1];
rz(-0.34698326) q[1];
sx q[1];
rz(0.68106236) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5446078) q[3];
sx q[3];
rz(-1.0506949) q[3];
sx q[3];
rz(2.6900684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3240933) q[2];
sx q[2];
rz(-1.9976511) q[2];
sx q[2];
rz(-0.49638805) q[2];
rz(3.056622) q[3];
sx q[3];
rz(-2.7111263) q[3];
sx q[3];
rz(2.5607204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519415) q[0];
sx q[0];
rz(-1.713151) q[0];
sx q[0];
rz(2.7570643) q[0];
rz(0.25228581) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(1.5581473) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88640942) q[0];
sx q[0];
rz(-0.46570581) q[0];
sx q[0];
rz(-3.1204771) q[0];
rz(-pi) q[1];
rz(-0.39790384) q[2];
sx q[2];
rz(-2.845484) q[2];
sx q[2];
rz(-2.1841911) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5024155) q[1];
sx q[1];
rz(-0.98160896) q[1];
sx q[1];
rz(-3.1261954) q[1];
rz(-2.4676314) q[3];
sx q[3];
rz(-2.1836851) q[3];
sx q[3];
rz(0.4919258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(-2.0015049) q[2];
rz(-1.7831066) q[3];
sx q[3];
rz(-1.8890817) q[3];
sx q[3];
rz(0.31806773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110905) q[0];
sx q[0];
rz(-1.6999812) q[0];
sx q[0];
rz(-0.39518133) q[0];
rz(-1.9197865) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(-2.2538259) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783274) q[0];
sx q[0];
rz(-2.3774231) q[0];
sx q[0];
rz(-2.7394638) q[0];
x q[1];
rz(1.0958065) q[2];
sx q[2];
rz(-1.57395) q[2];
sx q[2];
rz(-0.17660618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99988538) q[1];
sx q[1];
rz(-1.5453242) q[1];
sx q[1];
rz(1.9266121) q[1];
x q[2];
rz(0.10668175) q[3];
sx q[3];
rz(-1.8802207) q[3];
sx q[3];
rz(0.14313652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44783529) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(0.57175222) q[2];
rz(-1.5104431) q[3];
sx q[3];
rz(-1.7084728) q[3];
sx q[3];
rz(1.800764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1258662) q[0];
sx q[0];
rz(-1.8342352) q[0];
sx q[0];
rz(-2.9974708) q[0];
rz(-1.1852599) q[1];
sx q[1];
rz(-0.85694255) q[1];
sx q[1];
rz(-1.8640929) q[1];
rz(-2.754517) q[2];
sx q[2];
rz(-1.7995926) q[2];
sx q[2];
rz(2.8298701) q[2];
rz(0.41729144) q[3];
sx q[3];
rz(-2.4536316) q[3];
sx q[3];
rz(-1.3578547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
