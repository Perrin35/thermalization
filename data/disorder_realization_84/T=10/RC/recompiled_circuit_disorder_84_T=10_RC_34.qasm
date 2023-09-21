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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1140808) q[0];
sx q[0];
rz(-1.0923166) q[0];
sx q[0];
rz(1.6041243) q[0];
rz(-pi) q[1];
rz(-0.16205807) q[2];
sx q[2];
rz(-1.098212) q[2];
sx q[2];
rz(0.87640793) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.3996934) q[1];
sx q[1];
rz(-1.0642991) q[1];
sx q[1];
rz(0.472881) q[1];
x q[2];
rz(-1.4046895) q[3];
sx q[3];
rz(-1.7475834) q[3];
sx q[3];
rz(3.0151031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51241088) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(0.50981057) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(-2.0289039) q[0];
rz(-0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(1.8033093) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3801549) q[0];
sx q[0];
rz(-2.2614711) q[0];
sx q[0];
rz(-1.3209016) q[0];
rz(2.3035405) q[2];
sx q[2];
rz(-1.8997955) q[2];
sx q[2];
rz(2.7555562) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4900134) q[1];
sx q[1];
rz(-2.7499008) q[1];
sx q[1];
rz(-0.3117805) q[1];
rz(-2.5541833) q[3];
sx q[3];
rz(-1.7055159) q[3];
sx q[3];
rz(0.87232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.082836941) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.9632957) q[2];
rz(-0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96175471) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(-0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(2.8082074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51334914) q[0];
sx q[0];
rz(-1.7772563) q[0];
sx q[0];
rz(0.87135656) q[0];
x q[1];
rz(3.0079191) q[2];
sx q[2];
rz(-1.1867656) q[2];
sx q[2];
rz(-3.0093699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9874939) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(3.0497453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7942579) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(0.37924757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(-2.3392759) q[2];
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
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(2.5634649) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(-2.633458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80721891) q[0];
sx q[0];
rz(-0.59893543) q[0];
sx q[0];
rz(0.81143023) q[0];
x q[1];
rz(3.1122909) q[2];
sx q[2];
rz(-0.63819956) q[2];
sx q[2];
rz(-2.4820941) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9307738) q[1];
sx q[1];
rz(-0.79508077) q[1];
sx q[1];
rz(-2.7706183) q[1];
rz(2.1448137) q[3];
sx q[3];
rz(-2.814631) q[3];
sx q[3];
rz(0.60861482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4529139) q[2];
sx q[2];
rz(-2.8082509) q[2];
sx q[2];
rz(0.63344669) q[2];
rz(-0.59988919) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7552345) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(-2.9454943) q[0];
rz(1.8798937) q[1];
sx q[1];
rz(-2.321107) q[1];
sx q[1];
rz(2.0702147) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3065942) q[0];
sx q[0];
rz(-2.1822565) q[0];
sx q[0];
rz(0.5284662) q[0];
rz(-1.1159665) q[2];
sx q[2];
rz(-2.358846) q[2];
sx q[2];
rz(-2.6775529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(0.96697076) q[1];
x q[2];
rz(2.3134872) q[3];
sx q[3];
rz(-0.58907408) q[3];
sx q[3];
rz(2.7554054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.034996899) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(2.1595947) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-2.2976112) q[3];
sx q[3];
rz(1.8765607) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-0.53034267) q[0];
rz(-1.8416587) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(2.9249654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7055571) q[0];
sx q[0];
rz(-0.85058054) q[0];
sx q[0];
rz(0.35465045) q[0];
rz(-1.2712619) q[2];
sx q[2];
rz(-2.1795863) q[2];
sx q[2];
rz(1.8744206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0385598) q[1];
sx q[1];
rz(-1.9060887) q[1];
sx q[1];
rz(-3.1328939) q[1];
x q[2];
rz(-2.9495903) q[3];
sx q[3];
rz(-2.2322906) q[3];
sx q[3];
rz(1.1166752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65016046) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(-1.023863) q[2];
rz(0.8762382) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(-1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6525456) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(2.6126557) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(1.0891917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3952336) q[0];
sx q[0];
rz(-2.8440209) q[0];
sx q[0];
rz(1.5901106) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1955397) q[2];
sx q[2];
rz(-2.3240528) q[2];
sx q[2];
rz(-2.3600876) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84767064) q[1];
sx q[1];
rz(-1.5772181) q[1];
sx q[1];
rz(-0.033172219) q[1];
rz(-3.0508556) q[3];
sx q[3];
rz(-1.4416579) q[3];
sx q[3];
rz(1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0163394) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(1.4922173) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4847223) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(-0.86762506) q[0];
rz(3.0743657) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(0.19518383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282043) q[0];
sx q[0];
rz(-1.2767316) q[0];
sx q[0];
rz(1.78079) q[0];
rz(-pi) q[1];
rz(-0.40114258) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(-1.0912947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.118604) q[1];
sx q[1];
rz(-1.6123562) q[1];
sx q[1];
rz(1.1409345) q[1];
rz(2.1956452) q[3];
sx q[3];
rz(-0.56560707) q[3];
sx q[3];
rz(0.92029508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8937257) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-2.2360738) q[2];
rz(1.9838105) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(2.982443) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72682196) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(-2.1006405) q[0];
rz(0.078661593) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(0.35531607) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5268363) q[0];
sx q[0];
rz(-0.38072452) q[0];
sx q[0];
rz(-1.3044796) q[0];
rz(-pi) q[1];
rz(-2.9012868) q[2];
sx q[2];
rz(-1.086364) q[2];
sx q[2];
rz(-0.74953178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3231343) q[1];
sx q[1];
rz(-1.6272021) q[1];
sx q[1];
rz(1.2468546) q[1];
rz(-pi) q[2];
rz(0.60042419) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(-0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(1.5374373) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733317) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(-2.6955993) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(-2.840852) q[1];
rz(pi/2) q[2];
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
x q[1];
rz(2.5228595) q[2];
sx q[2];
rz(-0.64090568) q[2];
sx q[2];
rz(-0.59782019) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7584596) q[1];
sx q[1];
rz(-2.8499095) q[1];
sx q[1];
rz(3.0848857) q[1];
x q[2];
rz(-0.91046393) q[3];
sx q[3];
rz(-1.904389) q[3];
sx q[3];
rz(-3.0433082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5131502) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(2.002031) q[2];
rz(1.7906174) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(1.7383472) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(1.5691112) q[2];
sx q[2];
rz(-1.5148074) q[2];
sx q[2];
rz(-2.2956216) q[2];
rz(-2.8786447) q[3];
sx q[3];
rz(-0.96649747) q[3];
sx q[3];
rz(-1.4154712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
