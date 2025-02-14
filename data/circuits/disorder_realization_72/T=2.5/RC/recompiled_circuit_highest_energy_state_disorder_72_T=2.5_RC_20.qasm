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
rz(-2.3209651) q[0];
sx q[0];
rz(-2.4790915) q[0];
sx q[0];
rz(-2.878046) q[0];
rz(1.1072371) q[1];
sx q[1];
rz(4.4237408) q[1];
sx q[1];
rz(10.10034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5476108) q[0];
sx q[0];
rz(-1.2351961) q[0];
sx q[0];
rz(-0.006860126) q[0];
x q[1];
rz(-0.98904977) q[2];
sx q[2];
rz(-2.0768696) q[2];
sx q[2];
rz(-2.842234) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2517667) q[1];
sx q[1];
rz(-1.7210362) q[1];
sx q[1];
rz(-1.9637693) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1249902) q[3];
sx q[3];
rz(-2.5678291) q[3];
sx q[3];
rz(0.91780182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2519309) q[2];
sx q[2];
rz(-1.9193005) q[2];
sx q[2];
rz(-2.7737854) q[2];
rz(0.31428549) q[3];
sx q[3];
rz(-0.29044423) q[3];
sx q[3];
rz(0.17816003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3322068) q[0];
sx q[0];
rz(-3.0520913) q[0];
sx q[0];
rz(1.7820763) q[0];
rz(-1.2844194) q[1];
sx q[1];
rz(-1.1651243) q[1];
sx q[1];
rz(0.051609106) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2211014) q[0];
sx q[0];
rz(-1.3055542) q[0];
sx q[0];
rz(-0.18764253) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2003635) q[2];
sx q[2];
rz(-1.7631754) q[2];
sx q[2];
rz(0.34483257) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1866086) q[1];
sx q[1];
rz(-0.37395769) q[1];
sx q[1];
rz(-1.5723521) q[1];
rz(-2.6671131) q[3];
sx q[3];
rz(-1.9360376) q[3];
sx q[3];
rz(-2.7547835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.08903683) q[2];
sx q[2];
rz(-0.26036981) q[2];
sx q[2];
rz(0.54711771) q[2];
rz(-1.8735006) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(0.33856302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7608305) q[0];
sx q[0];
rz(-2.1050504) q[0];
sx q[0];
rz(-1.3408252) q[0];
rz(2.2876168) q[1];
sx q[1];
rz(-1.8546591) q[1];
sx q[1];
rz(2.8391848) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15203619) q[0];
sx q[0];
rz(-1.8468509) q[0];
sx q[0];
rz(2.2130475) q[0];
x q[1];
rz(2.6906784) q[2];
sx q[2];
rz(-1.779777) q[2];
sx q[2];
rz(-1.9318888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9309792) q[1];
sx q[1];
rz(-0.96081464) q[1];
sx q[1];
rz(-2.239834) q[1];
rz(2.0539927) q[3];
sx q[3];
rz(-0.54387605) q[3];
sx q[3];
rz(-2.9915006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93620187) q[2];
sx q[2];
rz(-1.8855636) q[2];
sx q[2];
rz(-2.5993247) q[2];
rz(-0.99644709) q[3];
sx q[3];
rz(-2.9139329) q[3];
sx q[3];
rz(-0.12797932) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19313136) q[0];
sx q[0];
rz(-1.4840935) q[0];
sx q[0];
rz(1.8782072) q[0];
rz(2.7073233) q[1];
sx q[1];
rz(-0.77249211) q[1];
sx q[1];
rz(-1.451452) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5562486) q[0];
sx q[0];
rz(-0.55216575) q[0];
sx q[0];
rz(-2.564173) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0201597) q[2];
sx q[2];
rz(-1.0673293) q[2];
sx q[2];
rz(1.891328) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5401104) q[1];
sx q[1];
rz(-0.14232351) q[1];
sx q[1];
rz(-0.92592923) q[1];
rz(-pi) q[2];
rz(2.1918815) q[3];
sx q[3];
rz(-1.1009645) q[3];
sx q[3];
rz(-2.7787152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0666535) q[2];
sx q[2];
rz(-0.94330072) q[2];
sx q[2];
rz(2.2554876) q[2];
rz(-3.1352299) q[3];
sx q[3];
rz(-1.3402091) q[3];
sx q[3];
rz(-2.0224915) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56115443) q[0];
sx q[0];
rz(-2.6851324) q[0];
sx q[0];
rz(1.2503257) q[0];
rz(2.9087032) q[1];
sx q[1];
rz(-0.75585514) q[1];
sx q[1];
rz(-3.0511391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6503432) q[0];
sx q[0];
rz(-1.8199662) q[0];
sx q[0];
rz(3.0155474) q[0];
x q[1];
rz(0.47201856) q[2];
sx q[2];
rz(-1.8119805) q[2];
sx q[2];
rz(2.3015353) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0098086987) q[1];
sx q[1];
rz(-2.4854641) q[1];
sx q[1];
rz(3.0387278) q[1];
rz(-pi) q[2];
rz(1.3227664) q[3];
sx q[3];
rz(-2.8042386) q[3];
sx q[3];
rz(-3.1243008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6470486) q[2];
sx q[2];
rz(-0.19379751) q[2];
sx q[2];
rz(-1.9336644) q[2];
rz(-2.2203994) q[3];
sx q[3];
rz(-0.84417206) q[3];
sx q[3];
rz(2.6827961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10814609) q[0];
sx q[0];
rz(-1.4220081) q[0];
sx q[0];
rz(-0.49149996) q[0];
rz(-0.72832251) q[1];
sx q[1];
rz(-1.1110543) q[1];
sx q[1];
rz(-2.947015) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74516836) q[0];
sx q[0];
rz(-1.0108507) q[0];
sx q[0];
rz(-0.26225175) q[0];
rz(-0.75677432) q[2];
sx q[2];
rz(-0.96454731) q[2];
sx q[2];
rz(-2.221125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42015612) q[1];
sx q[1];
rz(-2.2682574) q[1];
sx q[1];
rz(-0.8832133) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4846912) q[3];
sx q[3];
rz(-1.3852775) q[3];
sx q[3];
rz(2.3168062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3538094) q[2];
sx q[2];
rz(-1.0156735) q[2];
sx q[2];
rz(-0.0018399012) q[2];
rz(1.4932102) q[3];
sx q[3];
rz(-0.73072481) q[3];
sx q[3];
rz(-0.28436896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2581185) q[0];
sx q[0];
rz(-2.879877) q[0];
sx q[0];
rz(1.0650241) q[0];
rz(0.26271543) q[1];
sx q[1];
rz(-0.5717259) q[1];
sx q[1];
rz(-2.6782742) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1048038) q[0];
sx q[0];
rz(-2.1533009) q[0];
sx q[0];
rz(2.9534485) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1049785) q[2];
sx q[2];
rz(-2.2232703) q[2];
sx q[2];
rz(-0.62803113) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9874679) q[1];
sx q[1];
rz(-2.6374258) q[1];
sx q[1];
rz(-0.28480633) q[1];
x q[2];
rz(-0.024302146) q[3];
sx q[3];
rz(-1.4557456) q[3];
sx q[3];
rz(-2.3772511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.00039438417) q[2];
sx q[2];
rz(-1.1606777) q[2];
sx q[2];
rz(0.82994962) q[2];
rz(1.0478033) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(-0.1786264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.33775) q[0];
sx q[0];
rz(-2.1512845) q[0];
sx q[0];
rz(0.75482279) q[0];
rz(0.38914514) q[1];
sx q[1];
rz(-1.4578338) q[1];
sx q[1];
rz(0.39774242) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3869285) q[0];
sx q[0];
rz(-0.3873741) q[0];
sx q[0];
rz(-2.3088916) q[0];
rz(-pi) q[1];
rz(-1.7988458) q[2];
sx q[2];
rz(-2.1104276) q[2];
sx q[2];
rz(1.9613105) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.039835416) q[1];
sx q[1];
rz(-1.3164742) q[1];
sx q[1];
rz(-1.2891585) q[1];
rz(-pi) q[2];
rz(2.3944309) q[3];
sx q[3];
rz(-2.0748169) q[3];
sx q[3];
rz(0.51591831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3079181) q[2];
sx q[2];
rz(-2.8454744) q[2];
sx q[2];
rz(-0.31416848) q[2];
rz(1.2416154) q[3];
sx q[3];
rz(-1.5528468) q[3];
sx q[3];
rz(1.2392905) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.817713) q[0];
sx q[0];
rz(-0.68809026) q[0];
sx q[0];
rz(-2.896198) q[0];
rz(-1.2303526) q[1];
sx q[1];
rz(-3.0407258) q[1];
sx q[1];
rz(-2.6255677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300872) q[0];
sx q[0];
rz(-1.9507512) q[0];
sx q[0];
rz(1.7854693) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2788827) q[2];
sx q[2];
rz(-2.1893775) q[2];
sx q[2];
rz(1.0363621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83938767) q[1];
sx q[1];
rz(-0.47014648) q[1];
sx q[1];
rz(0.66580982) q[1];
rz(-pi) q[2];
rz(0.42428942) q[3];
sx q[3];
rz(-2.9919713) q[3];
sx q[3];
rz(-0.77009088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5838991) q[2];
sx q[2];
rz(-1.2217174) q[2];
sx q[2];
rz(-1.0830797) q[2];
rz(2.9205868) q[3];
sx q[3];
rz(-0.14328863) q[3];
sx q[3];
rz(-2.3451282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1204656) q[0];
sx q[0];
rz(-3.0607304) q[0];
sx q[0];
rz(3.0057111) q[0];
rz(-1.4295626) q[1];
sx q[1];
rz(-2.0202961) q[1];
sx q[1];
rz(2.8543465) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3843975) q[0];
sx q[0];
rz(-1.7585737) q[0];
sx q[0];
rz(1.8191992) q[0];
rz(-1.3332149) q[2];
sx q[2];
rz(-1.6529978) q[2];
sx q[2];
rz(1.0437708) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29746398) q[1];
sx q[1];
rz(-1.7904568) q[1];
sx q[1];
rz(2.541887) q[1];
rz(-pi) q[2];
rz(-0.29218896) q[3];
sx q[3];
rz(-2.1351519) q[3];
sx q[3];
rz(-0.60114229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91934943) q[2];
sx q[2];
rz(-1.2390077) q[2];
sx q[2];
rz(-0.83402056) q[2];
rz(2.833994) q[3];
sx q[3];
rz(-0.37720507) q[3];
sx q[3];
rz(0.26267499) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1766227) q[0];
sx q[0];
rz(-1.3055834) q[0];
sx q[0];
rz(2.1156043) q[0];
rz(2.7550244) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(1.4390527) q[2];
sx q[2];
rz(-1.9921039) q[2];
sx q[2];
rz(2.8767246) q[2];
rz(0.97053075) q[3];
sx q[3];
rz(-2.8875792) q[3];
sx q[3];
rz(-2.447788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
