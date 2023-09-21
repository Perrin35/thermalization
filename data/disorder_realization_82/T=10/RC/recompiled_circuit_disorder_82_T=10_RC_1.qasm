OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5582433) q[0];
sx q[0];
rz(-1.8177176) q[0];
sx q[0];
rz(-2.9980744) q[0];
x q[1];
rz(1.7945292) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(-1.7681233) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0919839) q[1];
sx q[1];
rz(-2.1503452) q[1];
sx q[1];
rz(2.0642573) q[1];
x q[2];
rz(0.7198556) q[3];
sx q[3];
rz(-0.30656439) q[3];
sx q[3];
rz(0.43678624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0063643) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(2.4025829) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(-3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(0.48746902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9092642) q[0];
sx q[0];
rz(-0.10350138) q[0];
sx q[0];
rz(2.2444025) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0080645) q[2];
sx q[2];
rz(-1.9568866) q[2];
sx q[2];
rz(-1.3646477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.060354787) q[1];
sx q[1];
rz(-1.758467) q[1];
sx q[1];
rz(2.9802122) q[1];
x q[2];
rz(1.8295733) q[3];
sx q[3];
rz(-2.8142455) q[3];
sx q[3];
rz(2.5395288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(1.0428492) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(-0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(2.9188459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73244625) q[0];
sx q[0];
rz(-1.961686) q[0];
sx q[0];
rz(-1.976165) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-0.97937102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0157156) q[1];
sx q[1];
rz(-1.5029969) q[1];
sx q[1];
rz(-2.5281639) q[1];
rz(-pi) q[2];
rz(0.95007105) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(1.0328968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(0.9202252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.675925) q[0];
sx q[0];
rz(-1.1629857) q[0];
sx q[0];
rz(2.2949335) q[0];
x q[1];
rz(2.9158981) q[2];
sx q[2];
rz(-2.7719471) q[2];
sx q[2];
rz(-0.29633488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0238266) q[1];
sx q[1];
rz(-0.95925602) q[1];
sx q[1];
rz(-2.141342) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3654877) q[3];
sx q[3];
rz(-1.9234386) q[3];
sx q[3];
rz(-2.7384788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-2.6884902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3081007) q[0];
sx q[0];
rz(-1.5335324) q[0];
sx q[0];
rz(1.4971855) q[0];
x q[1];
rz(-3.1066936) q[2];
sx q[2];
rz(-1.4279832) q[2];
sx q[2];
rz(-2.3226483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6140155) q[1];
sx q[1];
rz(-1.7909389) q[1];
sx q[1];
rz(0.084053587) q[1];
rz(-2.1702607) q[3];
sx q[3];
rz(-1.1820275) q[3];
sx q[3];
rz(2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-2.2912912) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033427514) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.5135117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244004) q[0];
sx q[0];
rz(-1.2478561) q[0];
sx q[0];
rz(0.49844235) q[0];
rz(-pi) q[1];
rz(1.1210576) q[2];
sx q[2];
rz(-1.5151086) q[2];
sx q[2];
rz(2.4186717) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.39735079) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(0.97126295) q[1];
x q[2];
rz(0.87168872) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(-2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.9980105) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.54388261) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163527) q[0];
sx q[0];
rz(-1.851474) q[0];
sx q[0];
rz(1.0513845) q[0];
x q[1];
rz(2.6830964) q[2];
sx q[2];
rz(-3.0756604) q[2];
sx q[2];
rz(-2.3022848) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.051441593) q[1];
sx q[1];
rz(-1.7045583) q[1];
sx q[1];
rz(0.83530463) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39077057) q[3];
sx q[3];
rz(-2.4739389) q[3];
sx q[3];
rz(-0.30927502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(1.7857893) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(1.0303248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.94901) q[0];
sx q[0];
rz(-2.2826676) q[0];
sx q[0];
rz(-2.683995) q[0];
x q[1];
rz(1.4872929) q[2];
sx q[2];
rz(-0.52831542) q[2];
sx q[2];
rz(-2.5022142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9955666) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(1.2703018) q[1];
rz(-pi) q[2];
rz(0.86352591) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(-1.2442148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(-2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-0.31627396) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-2.0064328) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.9599008) q[0];
sx q[0];
rz(-1.367021) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8912002) q[2];
sx q[2];
rz(-1.5431297) q[2];
sx q[2];
rz(0.49081341) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3865349) q[1];
sx q[1];
rz(-0.48663501) q[1];
sx q[1];
rz(1.6965894) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.023852392) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(2.9305487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(2.9454254) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553186) q[0];
sx q[0];
rz(-2.1325169) q[0];
sx q[0];
rz(2.4713211) q[0];
rz(0.050373366) q[2];
sx q[2];
rz(-1.7059687) q[2];
sx q[2];
rz(-2.2809575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4539459) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(0.83666283) q[1];
rz(-pi) q[2];
rz(-2.2979126) q[3];
sx q[3];
rz(-2.0920678) q[3];
sx q[3];
rz(2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(-3.0604559) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.1186196) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(1.6191471) q[2];
sx q[2];
rz(-1.3370677) q[2];
sx q[2];
rz(1.7481902) q[2];
rz(0.046394596) q[3];
sx q[3];
rz(-1.8002602) q[3];
sx q[3];
rz(0.044063448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];