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
rz(-2.2313843) q[0];
sx q[0];
rz(-2.3910523) q[0];
sx q[0];
rz(2.5810177) q[0];
rz(1.9380467) q[1];
sx q[1];
rz(-2.7653341) q[1];
sx q[1];
rz(-0.19317746) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3681741) q[0];
sx q[0];
rz(-2.3791398) q[0];
sx q[0];
rz(-2.4763251) q[0];
rz(-pi) q[1];
rz(1.4273663) q[2];
sx q[2];
rz(-0.71038112) q[2];
sx q[2];
rz(-1.3262891) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59483278) q[1];
sx q[1];
rz(-1.431219) q[1];
sx q[1];
rz(-2.2235893) q[1];
x q[2];
rz(-2.6704392) q[3];
sx q[3];
rz(-0.36590016) q[3];
sx q[3];
rz(0.19972502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91118497) q[2];
sx q[2];
rz(-2.513803) q[2];
sx q[2];
rz(0.5578624) q[2];
rz(-1.9134391) q[3];
sx q[3];
rz(-1.6817254) q[3];
sx q[3];
rz(-3.0465928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8519583) q[0];
sx q[0];
rz(-1.8486706) q[0];
sx q[0];
rz(1.5648382) q[0];
rz(1.9948888) q[1];
sx q[1];
rz(-1.7161938) q[1];
sx q[1];
rz(-2.1099405) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0609425) q[0];
sx q[0];
rz(-0.58851465) q[0];
sx q[0];
rz(1.1400162) q[0];
rz(-pi) q[1];
rz(-2.637219) q[2];
sx q[2];
rz(-2.6775355) q[2];
sx q[2];
rz(1.1936587) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90742753) q[1];
sx q[1];
rz(-1.9112559) q[1];
sx q[1];
rz(-2.2937066) q[1];
x q[2];
rz(-1.366721) q[3];
sx q[3];
rz(-2.0672247) q[3];
sx q[3];
rz(1.4816517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8138294) q[2];
sx q[2];
rz(-2.8507865) q[2];
sx q[2];
rz(1.4168868) q[2];
rz(-2.318577) q[3];
sx q[3];
rz(-1.5282642) q[3];
sx q[3];
rz(-0.64457646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7773892) q[0];
sx q[0];
rz(-2.0222029) q[0];
sx q[0];
rz(2.7893344) q[0];
rz(-0.38905713) q[1];
sx q[1];
rz(-1.16951) q[1];
sx q[1];
rz(-0.22638098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8181747) q[0];
sx q[0];
rz(-1.3238088) q[0];
sx q[0];
rz(-0.085121827) q[0];
x q[1];
rz(-0.34778123) q[2];
sx q[2];
rz(-1.0248597) q[2];
sx q[2];
rz(2.7394061) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4714053) q[1];
sx q[1];
rz(-0.67882631) q[1];
sx q[1];
rz(0.77969867) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6466633) q[3];
sx q[3];
rz(-1.1359204) q[3];
sx q[3];
rz(0.42805373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.160673) q[2];
sx q[2];
rz(-1.0272762) q[2];
sx q[2];
rz(-0.36160198) q[2];
rz(2.8412039) q[3];
sx q[3];
rz(-1.3114248) q[3];
sx q[3];
rz(0.96295199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.921973) q[0];
sx q[0];
rz(-2.0491056) q[0];
sx q[0];
rz(0.12119448) q[0];
rz(1.9425862) q[1];
sx q[1];
rz(-2.0620748) q[1];
sx q[1];
rz(1.8669063) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5124487) q[0];
sx q[0];
rz(-1.3487909) q[0];
sx q[0];
rz(-1.6187526) q[0];
rz(-pi) q[1];
x q[1];
rz(0.077353296) q[2];
sx q[2];
rz(-2.935545) q[2];
sx q[2];
rz(1.8872583) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.074305943) q[1];
sx q[1];
rz(-0.61903051) q[1];
sx q[1];
rz(2.7561042) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2043456) q[3];
sx q[3];
rz(-0.38903686) q[3];
sx q[3];
rz(-1.4938482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.049383) q[2];
sx q[2];
rz(-1.4468687) q[2];
sx q[2];
rz(1.764074) q[2];
rz(2.0130017) q[3];
sx q[3];
rz(-0.94213525) q[3];
sx q[3];
rz(-2.3142864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.2500797) q[0];
sx q[0];
rz(-2.2703607) q[0];
sx q[0];
rz(2.0569892) q[0];
rz(-2.4705823) q[1];
sx q[1];
rz(-1.6295461) q[1];
sx q[1];
rz(3.0674518) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0468354) q[0];
sx q[0];
rz(-1.0124542) q[0];
sx q[0];
rz(1.2461016) q[0];
x q[1];
rz(0.83168516) q[2];
sx q[2];
rz(-0.41809575) q[2];
sx q[2];
rz(1.3379793) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74846327) q[1];
sx q[1];
rz(-2.4824351) q[1];
sx q[1];
rz(1.1200957) q[1];
x q[2];
rz(0.13617326) q[3];
sx q[3];
rz(-1.2981997) q[3];
sx q[3];
rz(0.54823247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1159749) q[2];
sx q[2];
rz(-0.74883777) q[2];
sx q[2];
rz(2.1785114) q[2];
rz(1.3747619) q[3];
sx q[3];
rz(-1.129496) q[3];
sx q[3];
rz(0.53205427) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9328203) q[0];
sx q[0];
rz(-0.29173478) q[0];
sx q[0];
rz(-1.2579086) q[0];
rz(-2.3014297) q[1];
sx q[1];
rz(-1.9636619) q[1];
sx q[1];
rz(0.69200969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7888579) q[0];
sx q[0];
rz(-1.2458572) q[0];
sx q[0];
rz(-0.28577013) q[0];
rz(-pi) q[1];
rz(2.9553086) q[2];
sx q[2];
rz(-2.3669503) q[2];
sx q[2];
rz(-1.2828522) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5816725) q[1];
sx q[1];
rz(-0.85216235) q[1];
sx q[1];
rz(0.88345673) q[1];
rz(-1.0133842) q[3];
sx q[3];
rz(-1.0104826) q[3];
sx q[3];
rz(-2.2309365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.96364) q[2];
sx q[2];
rz(-1.0963564) q[2];
sx q[2];
rz(2.625551) q[2];
rz(-1.7959203) q[3];
sx q[3];
rz(-0.44107744) q[3];
sx q[3];
rz(1.0015944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26821414) q[0];
sx q[0];
rz(-0.28609797) q[0];
sx q[0];
rz(-2.1116665) q[0];
rz(0.65451199) q[1];
sx q[1];
rz(-1.3348568) q[1];
sx q[1];
rz(-2.9676504) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2692341) q[0];
sx q[0];
rz(-0.96251153) q[0];
sx q[0];
rz(-2.48003) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51474606) q[2];
sx q[2];
rz(-2.2651849) q[2];
sx q[2];
rz(1.0389164) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.084981266) q[1];
sx q[1];
rz(-0.95889303) q[1];
sx q[1];
rz(-2.8065744) q[1];
rz(-pi) q[2];
rz(-0.44486041) q[3];
sx q[3];
rz(-2.5813817) q[3];
sx q[3];
rz(2.1261393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30764636) q[2];
sx q[2];
rz(-2.6070194) q[2];
sx q[2];
rz(-1.8024811) q[2];
rz(0.80859679) q[3];
sx q[3];
rz(-1.9924889) q[3];
sx q[3];
rz(1.2861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22063743) q[0];
sx q[0];
rz(-2.0957102) q[0];
sx q[0];
rz(1.165423) q[0];
rz(2.9479345) q[1];
sx q[1];
rz(-0.8332738) q[1];
sx q[1];
rz(1.1509034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6100463) q[0];
sx q[0];
rz(-0.70554107) q[0];
sx q[0];
rz(0.39682589) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80457689) q[2];
sx q[2];
rz(-1.6940306) q[2];
sx q[2];
rz(-3.0446788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.833646) q[1];
sx q[1];
rz(-2.4599791) q[1];
sx q[1];
rz(2.2824085) q[1];
rz(0.99566136) q[3];
sx q[3];
rz(-2.1365949) q[3];
sx q[3];
rz(-0.034269661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9264441) q[2];
sx q[2];
rz(-2.1200659) q[2];
sx q[2];
rz(1.8683757) q[2];
rz(-2.6269954) q[3];
sx q[3];
rz(-2.8223346) q[3];
sx q[3];
rz(-2.4014421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.99106818) q[0];
sx q[0];
rz(-3.0758698) q[0];
sx q[0];
rz(-0.42732987) q[0];
rz(-1.7006251) q[1];
sx q[1];
rz(-1.7040375) q[1];
sx q[1];
rz(0.89868054) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81954992) q[0];
sx q[0];
rz(-1.2614095) q[0];
sx q[0];
rz(1.1197107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91941339) q[2];
sx q[2];
rz(-0.69910895) q[2];
sx q[2];
rz(1.0369911) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8391621) q[1];
sx q[1];
rz(-1.7453472) q[1];
sx q[1];
rz(1.8086834) q[1];
x q[2];
rz(-1.1077088) q[3];
sx q[3];
rz(-2.0631972) q[3];
sx q[3];
rz(-1.2295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23540363) q[2];
sx q[2];
rz(-1.0455422) q[2];
sx q[2];
rz(-1.7380627) q[2];
rz(-0.67052001) q[3];
sx q[3];
rz(-1.5454005) q[3];
sx q[3];
rz(-1.8286573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8197166) q[0];
sx q[0];
rz(-2.1724367) q[0];
sx q[0];
rz(-0.72625351) q[0];
rz(-2.4655474) q[1];
sx q[1];
rz(-1.36422) q[1];
sx q[1];
rz(-2.3647251) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60706106) q[0];
sx q[0];
rz(-0.93124572) q[0];
sx q[0];
rz(3.0846918) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6157584) q[2];
sx q[2];
rz(-1.4248214) q[2];
sx q[2];
rz(0.0073429664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90112359) q[1];
sx q[1];
rz(-1.2637584) q[1];
sx q[1];
rz(1.9528051) q[1];
rz(-pi) q[2];
rz(1.0105206) q[3];
sx q[3];
rz(-2.0891857) q[3];
sx q[3];
rz(-2.3862402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0016510222) q[2];
sx q[2];
rz(-0.17639128) q[2];
sx q[2];
rz(-1.4886935) q[2];
rz(-1.1267003) q[3];
sx q[3];
rz(-1.1688346) q[3];
sx q[3];
rz(-2.6039629) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61459944) q[0];
sx q[0];
rz(-2.499883) q[0];
sx q[0];
rz(-1.5747621) q[0];
rz(-0.78631403) q[1];
sx q[1];
rz(-0.17335261) q[1];
sx q[1];
rz(-0.39549624) q[1];
rz(-2.4509571) q[2];
sx q[2];
rz(-2.3403559) q[2];
sx q[2];
rz(1.5834119) q[2];
rz(1.8128822) q[3];
sx q[3];
rz(-2.3314706) q[3];
sx q[3];
rz(-1.7580845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
