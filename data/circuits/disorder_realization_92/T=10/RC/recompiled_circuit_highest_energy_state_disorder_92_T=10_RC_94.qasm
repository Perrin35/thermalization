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
rz(0.1306611) q[0];
sx q[0];
rz(4.9424439) q[0];
sx q[0];
rz(10.688936) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(-1.0882508) q[1];
sx q[1];
rz(0.01297125) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2269079) q[0];
sx q[0];
rz(-2.6400262) q[0];
sx q[0];
rz(-0.33059327) q[0];
rz(-pi) q[1];
rz(0.58836909) q[2];
sx q[2];
rz(-0.19489637) q[2];
sx q[2];
rz(2.4590059) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6889292) q[1];
sx q[1];
rz(-0.48349342) q[1];
sx q[1];
rz(-0.74093282) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8230328) q[3];
sx q[3];
rz(-2.5346642) q[3];
sx q[3];
rz(-1.3598809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26244792) q[2];
sx q[2];
rz(-1.6877354) q[2];
sx q[2];
rz(3.0461779) q[2];
rz(-0.48424193) q[3];
sx q[3];
rz(-2.3384194) q[3];
sx q[3];
rz(-1.6835083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478304) q[0];
sx q[0];
rz(-1.7050803) q[0];
sx q[0];
rz(2.4875212) q[0];
rz(1.3520799) q[1];
sx q[1];
rz(-0.65579954) q[1];
sx q[1];
rz(3.0316614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1760565) q[0];
sx q[0];
rz(-1.5829464) q[0];
sx q[0];
rz(-1.5506844) q[0];
rz(-pi) q[1];
rz(2.6835126) q[2];
sx q[2];
rz(-2.1501679) q[2];
sx q[2];
rz(-1.9922076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0866894) q[1];
sx q[1];
rz(-0.48799054) q[1];
sx q[1];
rz(-0.37505522) q[1];
rz(-2.7450493) q[3];
sx q[3];
rz(-1.8064883) q[3];
sx q[3];
rz(-1.9901708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.316651) q[2];
sx q[2];
rz(-1.2085088) q[2];
sx q[2];
rz(-1.6744772) q[2];
rz(-2.1064827) q[3];
sx q[3];
rz(-2.3605774) q[3];
sx q[3];
rz(-1.8234113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3159981) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(-2.3509534) q[0];
rz(-2.3979893) q[1];
sx q[1];
rz(-1.5433106) q[1];
sx q[1];
rz(0.57949439) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52344874) q[0];
sx q[0];
rz(-1.0001527) q[0];
sx q[0];
rz(2.4916925) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8786009) q[2];
sx q[2];
rz(-1.8557669) q[2];
sx q[2];
rz(-2.6927352) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.079536572) q[1];
sx q[1];
rz(-1.3128451) q[1];
sx q[1];
rz(3.1029812) q[1];
rz(-pi) q[2];
rz(2.0602588) q[3];
sx q[3];
rz(-2.0474985) q[3];
sx q[3];
rz(-0.34927327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4641331) q[2];
sx q[2];
rz(-2.3049998) q[2];
sx q[2];
rz(-0.38132384) q[2];
rz(-0.85121202) q[3];
sx q[3];
rz(-2.2150453) q[3];
sx q[3];
rz(2.4066431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6094991) q[0];
sx q[0];
rz(-0.61436009) q[0];
sx q[0];
rz(0.36439782) q[0];
rz(0.20026194) q[1];
sx q[1];
rz(-1.4118782) q[1];
sx q[1];
rz(3.0416378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17804201) q[0];
sx q[0];
rz(-0.077712312) q[0];
sx q[0];
rz(-1.8366209) q[0];
rz(1.3447138) q[2];
sx q[2];
rz(-2.7030309) q[2];
sx q[2];
rz(-1.7600702) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8980838) q[1];
sx q[1];
rz(-1.5484515) q[1];
sx q[1];
rz(-0.0925272) q[1];
x q[2];
rz(-2.5629879) q[3];
sx q[3];
rz(-2.1044995) q[3];
sx q[3];
rz(0.54909517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.072307) q[2];
sx q[2];
rz(-1.0682718) q[2];
sx q[2];
rz(3.0266673) q[2];
rz(1.140444) q[3];
sx q[3];
rz(-1.2830257) q[3];
sx q[3];
rz(0.97525245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9370148) q[0];
sx q[0];
rz(-2.1890722) q[0];
sx q[0];
rz(-0.17247795) q[0];
rz(-1.0109488) q[1];
sx q[1];
rz(-2.5806249) q[1];
sx q[1];
rz(0.15484658) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5342543) q[0];
sx q[0];
rz(-1.7887428) q[0];
sx q[0];
rz(-1.8740669) q[0];
x q[1];
rz(1.2482951) q[2];
sx q[2];
rz(-0.66993827) q[2];
sx q[2];
rz(-0.33316406) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1678007) q[1];
sx q[1];
rz(-1.1607045) q[1];
sx q[1];
rz(-1.4136259) q[1];
x q[2];
rz(-1.4590635) q[3];
sx q[3];
rz(-1.3936685) q[3];
sx q[3];
rz(-0.88645173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.054472063) q[2];
sx q[2];
rz(-0.27721578) q[2];
sx q[2];
rz(0.37230125) q[2];
rz(2.7436658) q[3];
sx q[3];
rz(-1.9979265) q[3];
sx q[3];
rz(-3.1411689) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.953124) q[0];
sx q[0];
rz(-2.5690881) q[0];
sx q[0];
rz(1.5931801) q[0];
rz(0.36987034) q[1];
sx q[1];
rz(-1.3984171) q[1];
sx q[1];
rz(2.5028548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8734147) q[0];
sx q[0];
rz(-1.333433) q[0];
sx q[0];
rz(-0.15837146) q[0];
rz(-2.0962786) q[2];
sx q[2];
rz(-1.7470226) q[2];
sx q[2];
rz(0.8187364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35013851) q[1];
sx q[1];
rz(-1.119507) q[1];
sx q[1];
rz(2.7103488) q[1];
rz(-pi) q[2];
rz(0.79213284) q[3];
sx q[3];
rz(-1.9895305) q[3];
sx q[3];
rz(-1.2611539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0686191) q[2];
sx q[2];
rz(-1.0516473) q[2];
sx q[2];
rz(-2.8886786) q[2];
rz(-2.2933293) q[3];
sx q[3];
rz(-1.4784644) q[3];
sx q[3];
rz(-3.0566791) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10097583) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(0.14895359) q[0];
rz(1.8303998) q[1];
sx q[1];
rz(-1.7313749) q[1];
sx q[1];
rz(-0.054100903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.427092) q[0];
sx q[0];
rz(-0.50760554) q[0];
sx q[0];
rz(-0.18113329) q[0];
rz(-0.80246822) q[2];
sx q[2];
rz(-1.6020498) q[2];
sx q[2];
rz(2.5078338) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5230519) q[1];
sx q[1];
rz(-0.7071916) q[1];
sx q[1];
rz(-0.16945736) q[1];
rz(-pi) q[2];
rz(1.6160433) q[3];
sx q[3];
rz(-2.6658969) q[3];
sx q[3];
rz(0.078381009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80186239) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(0.96699634) q[2];
rz(2.4407834) q[3];
sx q[3];
rz(-0.98027027) q[3];
sx q[3];
rz(-0.35082671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12857777) q[0];
sx q[0];
rz(-1.1739434) q[0];
sx q[0];
rz(-1.5933734) q[0];
rz(2.7633527) q[1];
sx q[1];
rz(-2.389237) q[1];
sx q[1];
rz(0.38984782) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10791099) q[0];
sx q[0];
rz(-1.0709239) q[0];
sx q[0];
rz(1.7938016) q[0];
rz(-pi) q[1];
rz(3.1249525) q[2];
sx q[2];
rz(-1.8748054) q[2];
sx q[2];
rz(-1.9250411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2521542) q[1];
sx q[1];
rz(-2.5233626) q[1];
sx q[1];
rz(-2.5095224) q[1];
rz(-0.71487074) q[3];
sx q[3];
rz(-0.33239543) q[3];
sx q[3];
rz(0.22163342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0491911) q[2];
sx q[2];
rz(-2.0202049) q[2];
sx q[2];
rz(1.258705) q[2];
rz(-2.8973268) q[3];
sx q[3];
rz(-1.3552856) q[3];
sx q[3];
rz(2.145483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.32023892) q[0];
sx q[0];
rz(-1.9122253) q[0];
sx q[0];
rz(1.7852596) q[0];
rz(-2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(0.43620268) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1937418) q[0];
sx q[0];
rz(-1.1202888) q[0];
sx q[0];
rz(-1.9175025) q[0];
rz(1.1567843) q[2];
sx q[2];
rz(-1.2591259) q[2];
sx q[2];
rz(3.0729635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0862941) q[1];
sx q[1];
rz(-2.3196967) q[1];
sx q[1];
rz(0.31676745) q[1];
rz(2.7049191) q[3];
sx q[3];
rz(-1.6323576) q[3];
sx q[3];
rz(1.3832645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0210375) q[2];
sx q[2];
rz(-1.1641116) q[2];
sx q[2];
rz(0.54350054) q[2];
rz(0.049526878) q[3];
sx q[3];
rz(-1.6560358) q[3];
sx q[3];
rz(1.653695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307584) q[0];
sx q[0];
rz(-0.13228358) q[0];
sx q[0];
rz(-3.0623867) q[0];
rz(0.40009701) q[1];
sx q[1];
rz(-0.85405093) q[1];
sx q[1];
rz(-2.1150463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.067018) q[0];
sx q[0];
rz(-1.3368784) q[0];
sx q[0];
rz(2.3090655) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8978149) q[2];
sx q[2];
rz(-0.37636435) q[2];
sx q[2];
rz(1.3805022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77791926) q[1];
sx q[1];
rz(-1.8936186) q[1];
sx q[1];
rz(1.670091) q[1];
rz(1.4913358) q[3];
sx q[3];
rz(-1.9206408) q[3];
sx q[3];
rz(-2.3235278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2058699) q[2];
sx q[2];
rz(-1.8368072) q[2];
sx q[2];
rz(-1.4523466) q[2];
rz(2.3116889) q[3];
sx q[3];
rz(-0.87703505) q[3];
sx q[3];
rz(-0.95480603) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6185388) q[0];
sx q[0];
rz(-1.9970311) q[0];
sx q[0];
rz(1.507623) q[0];
rz(1.2670831) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(-2.3054068) q[2];
sx q[2];
rz(-0.7291353) q[2];
sx q[2];
rz(1.8020204) q[2];
rz(-0.51658705) q[3];
sx q[3];
rz(-2.5174375) q[3];
sx q[3];
rz(0.843912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
