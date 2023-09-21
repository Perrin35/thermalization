OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(-2.1043632) q[0];
sx q[0];
rz(-0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1188388) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(1.3214146) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6718164) q[2];
sx q[2];
rz(-0.46576408) q[2];
sx q[2];
rz(2.2848406) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76525926) q[1];
sx q[1];
rz(-1.9783124) q[1];
sx q[1];
rz(0.63912649) q[1];
rz(-pi) q[2];
rz(-2.9079307) q[3];
sx q[3];
rz(-1.3705001) q[3];
sx q[3];
rz(1.3113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(0.65650666) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-0.48746902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66747626) q[0];
sx q[0];
rz(-1.6352909) q[0];
sx q[0];
rz(-1.4897896) q[0];
rz(-pi) q[1];
rz(-2.3357453) q[2];
sx q[2];
rz(-2.5666551) q[2];
sx q[2];
rz(0.88428674) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5408052) q[1];
sx q[1];
rz(-1.729319) q[1];
sx q[1];
rz(1.7608789) q[1];
rz(1.8879714) q[3];
sx q[3];
rz(-1.653169) q[3];
sx q[3];
rz(1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(-0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(2.9512067) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(-0.22274676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091464) q[0];
sx q[0];
rz(-1.1799066) q[0];
sx q[0];
rz(-1.1654277) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-2.1622216) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.682444) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(-0.1174121) q[1];
rz(-1.907903) q[3];
sx q[3];
rz(-1.8030231) q[3];
sx q[3];
rz(3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(-2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(-2.1779493) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(0.9202252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3738149) q[0];
sx q[0];
rz(-0.91705634) q[0];
sx q[0];
rz(2.6184404) q[0];
x q[1];
rz(1.6572861) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(2.6038225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0238266) q[1];
sx q[1];
rz(-0.95925602) q[1];
sx q[1];
rz(1.0002506) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50587378) q[3];
sx q[3];
rz(-0.40588356) q[3];
sx q[3];
rz(-0.94569262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9437287) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.7368332) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(2.6884902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720854) q[0];
sx q[0];
rz(-0.082490248) q[0];
sx q[0];
rz(-2.0399658) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7136953) q[2];
sx q[2];
rz(-1.60534) q[2];
sx q[2];
rz(-0.75682109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.52757712) q[1];
sx q[1];
rz(-1.3506538) q[1];
sx q[1];
rz(0.084053587) q[1];
x q[2];
rz(-2.681053) q[3];
sx q[3];
rz(-2.1201049) q[3];
sx q[3];
rz(2.3159546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0129464) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(0.72367469) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(0.85030142) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(1.628081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3753525) q[0];
sx q[0];
rz(-1.1002812) q[0];
sx q[0];
rz(-1.9348295) q[0];
rz(-pi) q[1];
rz(-1.4432625) q[2];
sx q[2];
rz(-2.6886534) q[2];
sx q[2];
rz(-2.4085101) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5634798) q[1];
sx q[1];
rz(-2.0260749) q[1];
sx q[1];
rz(-0.77225765) q[1];
rz(-pi) q[2];
rz(2.7558277) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(-0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.05904077) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(2.7440199) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.9546753) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2448954) q[0];
sx q[0];
rz(-0.58422409) q[0];
sx q[0];
rz(-2.0969735) q[0];
x q[1];
rz(-3.0824532) q[2];
sx q[2];
rz(-1.5999609) q[2];
sx q[2];
rz(-1.9524706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.051441593) q[1];
sx q[1];
rz(-1.7045583) q[1];
sx q[1];
rz(-2.306288) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7508221) q[3];
sx q[3];
rz(-2.4739389) q[3];
sx q[3];
rz(-2.8323176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.879803) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(1.5073744) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-2.2858455) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5949769) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(-2.0440408) q[0];
rz(3.0929504) q[2];
sx q[2];
rz(-2.0970793) q[2];
sx q[2];
rz(-2.4056048) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9955666) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(1.8712908) q[1];
x q[2];
rz(0.97708948) q[3];
sx q[3];
rz(-1.1247375) q[3];
sx q[3];
rz(-2.2462728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29256233) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(0.31627396) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(2.0064328) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0907744) q[0];
sx q[0];
rz(-2.7047815) q[0];
sx q[0];
rz(-0.45849053) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8912002) q[2];
sx q[2];
rz(-1.5431297) q[2];
sx q[2];
rz(-2.6507792) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75505776) q[1];
sx q[1];
rz(-2.6549576) q[1];
sx q[1];
rz(1.6965894) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1177403) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(-0.21104392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90821663) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(-2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88403945) q[0];
sx q[0];
rz(-2.1243874) q[0];
sx q[0];
rz(-2.2474225) q[0];
rz(-pi) q[1];
rz(-1.4354544) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(-0.70336715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6876467) q[1];
sx q[1];
rz(-1.5861662) q[1];
sx q[1];
rz(-2.3049298) q[1];
x q[2];
rz(0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(-2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3907884) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(0.081136726) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(-1.8267869) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(2.9076004) q[2];
sx q[2];
rz(-1.5237612) q[2];
sx q[2];
rz(0.16618726) q[2];
rz(-1.3410939) q[3];
sx q[3];
rz(-1.5256186) q[3];
sx q[3];
rz(-1.516173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];