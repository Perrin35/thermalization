OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8755655) q[0];
sx q[0];
rz(-2.6063483) q[0];
sx q[0];
rz(2.3875561) q[0];
rz(0.79020483) q[1];
sx q[1];
rz(-1.2269998) q[1];
sx q[1];
rz(-1.1608646) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95522674) q[0];
sx q[0];
rz(-0.47954924) q[0];
sx q[0];
rz(3.0774374) q[0];
rz(-pi) q[1];
rz(2.0487469) q[2];
sx q[2];
rz(-1.4266326) q[2];
sx q[2];
rz(2.521487) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7418993) q[1];
sx q[1];
rz(-2.0772935) q[1];
sx q[1];
rz(-2.6687117) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17920223) q[3];
sx q[3];
rz(-1.7342907) q[3];
sx q[3];
rz(-1.6678099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51241088) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(0.50981057) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(-1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78012413) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(1.1126888) q[0];
rz(2.9878222) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(1.3382834) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7614377) q[0];
sx q[0];
rz(-2.2614711) q[0];
sx q[0];
rz(1.3209016) q[0];
x q[1];
rz(-0.43054994) q[2];
sx q[2];
rz(-2.2562648) q[2];
sx q[2];
rz(-1.6738883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3160623) q[1];
sx q[1];
rz(-1.1989374) q[1];
sx q[1];
rz(1.4447681) q[1];
x q[2];
rz(-1.7322147) q[3];
sx q[3];
rz(-0.98940778) q[3];
sx q[3];
rz(0.60928173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.082836941) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(1.9632957) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798379) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(0.33338526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282368) q[0];
sx q[0];
rz(-0.88909273) q[0];
sx q[0];
rz(0.26716726) q[0];
rz(-pi) q[1];
rz(-1.9579499) q[2];
sx q[2];
rz(-1.4469115) q[2];
sx q[2];
rz(-1.7533592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15409878) q[1];
sx q[1];
rz(-2.0425659) q[1];
sx q[1];
rz(3.0497453) q[1];
x q[2];
rz(-2.5081627) q[3];
sx q[3];
rz(-1.7519577) q[3];
sx q[3];
rz(1.0599979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(0.80231673) q[2];
rz(-0.24117593) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(0.10087092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872831) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(2.5774082) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(-2.633458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6619686) q[0];
sx q[0];
rz(-1.9920252) q[0];
sx q[0];
rz(-0.43930102) q[0];
rz(-pi) q[1];
rz(3.1122909) q[2];
sx q[2];
rz(-0.63819956) q[2];
sx q[2];
rz(-2.4820941) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21081884) q[1];
sx q[1];
rz(-0.79508077) q[1];
sx q[1];
rz(2.7706183) q[1];
x q[2];
rz(-0.18211256) q[3];
sx q[3];
rz(-1.2977227) q[3];
sx q[3];
rz(1.9338716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(-2.5417035) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863581) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(-2.9454943) q[0];
rz(-1.8798937) q[1];
sx q[1];
rz(-2.321107) q[1];
sx q[1];
rz(1.0713779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51273726) q[0];
sx q[0];
rz(-2.3561986) q[0];
sx q[0];
rz(-0.94731713) q[0];
rz(-pi) q[1];
rz(-1.1159665) q[2];
sx q[2];
rz(-2.358846) q[2];
sx q[2];
rz(0.46403971) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86712671) q[1];
sx q[1];
rz(-1.9879838) q[1];
sx q[1];
rz(2.8309114) q[1];
rz(2.7171633) q[3];
sx q[3];
rz(-1.1491346) q[3];
sx q[3];
rz(1.9205586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(0.98199797) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035325) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(2.61125) q[0];
rz(1.2999339) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(-0.21662724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43603555) q[0];
sx q[0];
rz(-0.85058054) q[0];
sx q[0];
rz(-2.7869422) q[0];
rz(-0.40041133) q[2];
sx q[2];
rz(-2.4715804) q[2];
sx q[2];
rz(2.3695721) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4649012) q[1];
sx q[1];
rz(-1.5790107) q[1];
sx q[1];
rz(1.2354922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3304177) q[3];
sx q[3];
rz(-2.4568395) q[3];
sx q[3];
rz(1.4231589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65016046) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(-1.023863) q[2];
rz(-2.2653545) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(-1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6525456) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(-0.52893692) q[0];
rz(-1.6128929) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(-1.0891917) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.335621) q[0];
sx q[0];
rz(-1.5764589) q[0];
sx q[0];
rz(1.2732768) q[0];
rz(-pi) q[1];
rz(0.85765526) q[2];
sx q[2];
rz(-2.0115888) q[2];
sx q[2];
rz(-2.8105274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84767064) q[1];
sx q[1];
rz(-1.5772181) q[1];
sx q[1];
rz(-3.1084204) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0508556) q[3];
sx q[3];
rz(-1.4416579) q[3];
sx q[3];
rz(1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12525325) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(-1.2255229) q[2];
rz(1.4922173) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(-2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4847223) q[0];
sx q[0];
rz(-0.06198922) q[0];
sx q[0];
rz(-0.86762506) q[0];
rz(-3.0743657) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(-2.9464088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.09571) q[0];
sx q[0];
rz(-1.7716496) q[0];
sx q[0];
rz(-2.8413089) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7404501) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(1.0912947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.118604) q[1];
sx q[1];
rz(-1.5292364) q[1];
sx q[1];
rz(2.0006581) q[1];
rz(0.35555367) q[3];
sx q[3];
rz(-2.0204633) q[3];
sx q[3];
rz(0.21330968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.247867) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(2.2360738) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(-2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4147707) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(1.0409521) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-2.9610596) q[1];
sx q[1];
rz(-2.7862766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040877) q[0];
sx q[0];
rz(-1.4728439) q[0];
sx q[0];
rz(1.9393001) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9955194) q[2];
sx q[2];
rz(-2.6051084) q[2];
sx q[2];
rz(1.9076965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73341093) q[1];
sx q[1];
rz(-1.8942041) q[1];
sx q[1];
rz(3.0820993) q[1];
rz(0.60042419) q[3];
sx q[3];
rz(-1.416559) q[3];
sx q[3];
rz(-2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4650402) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733317) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(-0.4459933) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(-0.30074063) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7556831) q[0];
sx q[0];
rz(-1.798238) q[0];
sx q[0];
rz(1.0230416) q[0];
x q[1];
rz(-1.9791331) q[2];
sx q[2];
rz(-2.0795341) q[2];
sx q[2];
rz(-3.0131154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2749394) q[1];
sx q[1];
rz(-1.5544974) q[1];
sx q[1];
rz(-0.29124041) q[1];
x q[2];
rz(-2.2311287) q[3];
sx q[3];
rz(-1.2372036) q[3];
sx q[3];
rz(-3.0433082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62844244) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(1.7906174) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6931077) q[0];
sx q[0];
rz(-1.2379452) q[0];
sx q[0];
rz(-2.2647279) q[0];
rz(-1.7383472) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(3.0856037) q[2];
sx q[2];
rz(-1.5724788) q[2];
sx q[2];
rz(2.4168617) q[2];
rz(0.26294796) q[3];
sx q[3];
rz(-0.96649747) q[3];
sx q[3];
rz(-1.4154712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
