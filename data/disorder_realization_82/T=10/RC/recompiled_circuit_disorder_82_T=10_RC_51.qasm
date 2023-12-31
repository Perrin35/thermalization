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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1188388) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(-1.820178) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7202397) q[2];
sx q[2];
rz(-1.3660649) q[2];
sx q[2];
rz(-2.8533964) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3763334) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(2.5024662) q[1];
x q[2];
rz(1.7765331) q[3];
sx q[3];
rz(-1.3418901) q[3];
sx q[3];
rz(-2.8347901) q[3];
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
rz(-2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-2.6541236) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23232846) q[0];
sx q[0];
rz(-0.10350138) q[0];
sx q[0];
rz(-0.89719015) q[0];
x q[1];
rz(2.3357453) q[2];
sx q[2];
rz(-2.5666551) q[2];
sx q[2];
rz(-0.88428674) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4844955) q[1];
sx q[1];
rz(-0.24689455) q[1];
sx q[1];
rz(2.2730278) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2536212) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(1.9272643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-2.7462192) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(2.9188459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4091464) q[0];
sx q[0];
rz(-1.1799066) q[0];
sx q[0];
rz(1.1654277) q[0];
rz(-pi) q[1];
rz(0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-2.1622216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60274117) q[1];
sx q[1];
rz(-0.95898421) q[1];
sx q[1];
rz(-1.6536504) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1915216) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(-1.4032723) q[0];
rz(1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-2.2213675) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677778) q[0];
sx q[0];
rz(-0.91705634) q[0];
sx q[0];
rz(-2.6184404) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6572861) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(2.6038225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4184119) q[1];
sx q[1];
rz(-0.81058093) q[1];
sx q[1];
rz(0.65631887) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7761049) q[3];
sx q[3];
rz(-1.9234386) q[3];
sx q[3];
rz(0.40311381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.2801923) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(-1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.7368332) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(2.6884902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720854) q[0];
sx q[0];
rz(-0.082490248) q[0];
sx q[0];
rz(-2.0399658) q[0];
rz(-1.7136953) q[2];
sx q[2];
rz(-1.5362527) q[2];
sx q[2];
rz(-0.75682109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52757712) q[1];
sx q[1];
rz(-1.7909389) q[1];
sx q[1];
rz(-3.0575391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.971332) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(-0.72367469) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(2.2912912) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1081651) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(-0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.628081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244004) q[0];
sx q[0];
rz(-1.8937366) q[0];
sx q[0];
rz(-0.49844235) q[0];
x q[1];
rz(-0.061821826) q[2];
sx q[2];
rz(-1.1218058) q[2];
sx q[2];
rz(-0.87473727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43135333) q[1];
sx q[1];
rz(-0.87190404) q[1];
sx q[1];
rz(2.5297574) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87168872) q[3];
sx q[3];
rz(-0.5629493) q[3];
sx q[3];
rz(-2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(1.9980105) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(0.17091621) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(-1.1869173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12524) q[0];
sx q[0];
rz(-1.2901187) q[0];
sx q[0];
rz(2.0902082) q[0];
rz(0.059139472) q[2];
sx q[2];
rz(-1.5999609) q[2];
sx q[2];
rz(1.189122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6656875) q[1];
sx q[1];
rz(-2.396282) q[1];
sx q[1];
rz(-1.3728632) q[1];
x q[2];
rz(-2.5116634) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.7857893) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(-0.021727173) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-1.0303248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54661575) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(-1.0975518) q[0];
rz(1.4872929) q[2];
sx q[2];
rz(-2.6132772) q[2];
sx q[2];
rz(2.5022142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6701263) q[1];
sx q[1];
rz(-1.2738436) q[1];
sx q[1];
rz(2.9832277) q[1];
rz(-0.86352591) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(-1.8973779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29256233) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(-0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(0.58644811) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(-0.31627396) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(1.1351599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.1816918) q[0];
sx q[0];
rz(-1.7745716) q[0];
rz(-pi) q[1];
rz(1.8912002) q[2];
sx q[2];
rz(-1.5431297) q[2];
sx q[2];
rz(2.6507792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3865349) q[1];
sx q[1];
rz(-0.48663501) q[1];
sx q[1];
rz(-1.4450032) q[1];
rz(-0.023852392) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(2.9305487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(-0.28963447) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90821663) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10712121) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(0.79191533) q[0];
rz(-0.050373366) q[2];
sx q[2];
rz(-1.435624) q[2];
sx q[2];
rz(0.86063517) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.13388053) q[1];
sx q[1];
rz(-2.4073283) q[1];
sx q[1];
rz(-1.5937362) q[1];
rz(-pi) q[2];
rz(-2.2979126) q[3];
sx q[3];
rz(-2.0920678) q[3];
sx q[3];
rz(-0.47714864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3907884) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186196) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(2.9076004) q[2];
sx q[2];
rz(-1.5237612) q[2];
sx q[2];
rz(0.16618726) q[2];
rz(-1.374791) q[3];
sx q[3];
rz(-0.23402611) q[3];
sx q[3];
rz(0.2454161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
