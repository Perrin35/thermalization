OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(-0.53524435) q[0];
sx q[0];
rz(-2.3875561) q[0];
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(-1.9807281) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1863659) q[0];
sx q[0];
rz(-2.6620434) q[0];
sx q[0];
rz(0.064155302) q[0];
rz(1.0928458) q[2];
sx q[2];
rz(-1.4266326) q[2];
sx q[2];
rz(0.62010566) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7418993) q[1];
sx q[1];
rz(-2.0772935) q[1];
sx q[1];
rz(0.472881) q[1];
x q[2];
rz(-0.74673064) q[3];
sx q[3];
rz(-0.24198469) q[3];
sx q[3];
rz(2.5065968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51241088) q[2];
sx q[2];
rz(-1.2966195) q[2];
sx q[2];
rz(-2.4751002) q[2];
rz(2.6317821) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(-1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78012413) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(1.1126888) q[0];
rz(-2.9878222) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(-1.8033093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7897658) q[0];
sx q[0];
rz(-1.7625945) q[0];
sx q[0];
rz(0.70621323) q[0];
rz(-pi) q[1];
rz(2.0427225) q[2];
sx q[2];
rz(-2.3510691) q[2];
sx q[2];
rz(-2.3015442) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9328623) q[1];
sx q[1];
rz(-1.4534229) q[1];
sx q[1];
rz(-0.37456234) q[1];
x q[2];
rz(-0.23985858) q[3];
sx q[3];
rz(-0.60088241) q[3];
sx q[3];
rz(2.2440653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.082836941) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(1.9632957) q[2];
rz(2.1792049) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798379) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(-1.5270365) q[0];
rz(2.4987192) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(-2.8082074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282368) q[0];
sx q[0];
rz(-2.2524999) q[0];
sx q[0];
rz(-0.26716726) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8893759) q[2];
sx q[2];
rz(-0.40553667) q[2];
sx q[2];
rz(0.47682724) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6830605) q[1];
sx q[1];
rz(-1.6525869) q[1];
sx q[1];
rz(-1.0973147) q[1];
rz(-0.63342996) q[3];
sx q[3];
rz(-1.7519577) q[3];
sx q[3];
rz(-1.0599979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(2.3392759) q[2];
rz(0.24117593) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-2.5774082) q[0];
rz(-2.5634649) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(0.50813466) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.479624) q[0];
sx q[0];
rz(-1.1495674) q[0];
sx q[0];
rz(-0.43930102) q[0];
x q[1];
rz(1.5925243) q[2];
sx q[2];
rz(-2.2086775) q[2];
sx q[2];
rz(2.4456172) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9307738) q[1];
sx q[1];
rz(-0.79508077) q[1];
sx q[1];
rz(-0.37097431) q[1];
rz(0.18211256) q[3];
sx q[3];
rz(-1.84387) q[3];
sx q[3];
rz(-1.207721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(0.59988919) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863581) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(2.9454943) q[0];
rz(-1.8798937) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-1.0713779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51273726) q[0];
sx q[0];
rz(-0.78539408) q[0];
sx q[0];
rz(-2.1942755) q[0];
x q[1];
rz(0.84153701) q[2];
sx q[2];
rz(-1.255799) q[2];
sx q[2];
rz(2.3685761) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5372249) q[1];
sx q[1];
rz(-0.51465263) q[1];
sx q[1];
rz(0.96697076) q[1];
rz(-pi) q[2];
rz(-0.82810546) q[3];
sx q[3];
rz(-0.58907408) q[3];
sx q[3];
rz(-0.38618726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(2.1595947) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-0.53034267) q[0];
rz(1.8416587) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(-0.21662724) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.066031) q[0];
sx q[0];
rz(-2.3530585) q[0];
sx q[0];
rz(-1.1939474) q[0];
rz(-pi) q[1];
rz(-1.8703307) q[2];
sx q[2];
rz(-2.1795863) q[2];
sx q[2];
rz(1.2671721) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10303282) q[1];
sx q[1];
rz(-1.235504) q[1];
sx q[1];
rz(-3.1328939) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3304177) q[3];
sx q[3];
rz(-0.68475311) q[3];
sx q[3];
rz(1.4231589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4914322) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(-2.2653545) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4890471) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(2.6126557) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(-1.0891917) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3952336) q[0];
sx q[0];
rz(-0.2975718) q[0];
sx q[0];
rz(1.5901106) q[0];
x q[1];
rz(2.1955397) q[2];
sx q[2];
rz(-0.81753987) q[2];
sx q[2];
rz(2.3600876) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91428077) q[1];
sx q[1];
rz(-0.033787878) q[1];
sx q[1];
rz(0.19126161) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7004622) q[3];
sx q[3];
rz(-1.4808169) q[3];
sx q[3];
rz(-2.9699096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(-1.9160697) q[2];
rz(1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.65687031) q[0];
sx q[0];
rz(-0.06198922) q[0];
sx q[0];
rz(0.86762506) q[0];
rz(-3.0743657) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(0.19518383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0458826) q[0];
sx q[0];
rz(-1.3699431) q[0];
sx q[0];
rz(0.30028371) q[0];
rz(-0.40114258) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(2.0502979) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.118604) q[1];
sx q[1];
rz(-1.6123562) q[1];
sx q[1];
rz(1.1409345) q[1];
rz(-pi) q[2];
rz(0.35555367) q[3];
sx q[3];
rz(-1.1211294) q[3];
sx q[3];
rz(-0.21330968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8937257) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(0.90551886) q[2];
rz(1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(-2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72682196) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(2.1006405) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(2.7862766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040877) q[0];
sx q[0];
rz(-1.4728439) q[0];
sx q[0];
rz(1.9393001) q[0];
rz(1.0742498) q[2];
sx q[2];
rz(-1.3585919) q[2];
sx q[2];
rz(0.70763904) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91868672) q[1];
sx q[1];
rz(-2.8129473) q[1];
sx q[1];
rz(-1.3952284) q[1];
rz(-pi) q[2];
rz(0.60042419) q[3];
sx q[3];
rz(-1.416559) q[3];
sx q[3];
rz(-2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(-1.4036277) q[2];
rz(0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86826098) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(-0.4459933) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(-2.840852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8200127) q[0];
sx q[0];
rz(-2.1029148) q[0];
sx q[0];
rz(2.8768455) q[0];
rz(-pi) q[1];
rz(1.1624596) q[2];
sx q[2];
rz(-2.0795341) q[2];
sx q[2];
rz(-3.0131154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69925752) q[1];
sx q[1];
rz(-1.8619969) q[1];
sx q[1];
rz(1.553781) q[1];
x q[2];
rz(1.0565287) q[3];
sx q[3];
rz(-2.4132055) q[3];
sx q[3];
rz(1.8715093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.62844244) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(2.002031) q[2];
rz(1.3509753) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(1.1736419) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6931077) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(1.7383472) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(-0.055988978) q[2];
sx q[2];
rz(-1.5724788) q[2];
sx q[2];
rz(2.4168617) q[2];
rz(1.9308405) q[3];
sx q[3];
rz(-2.4891709) q[3];
sx q[3];
rz(-1.8579033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];