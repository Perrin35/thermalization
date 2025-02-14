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
rz(0.91520619) q[0];
sx q[0];
rz(-0.45867607) q[0];
sx q[0];
rz(0.31804481) q[0];
rz(1.3660499) q[1];
sx q[1];
rz(-2.4429758) q[1];
sx q[1];
rz(0.06279343) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013115766) q[0];
sx q[0];
rz(-2.4139691) q[0];
sx q[0];
rz(0.81650556) q[0];
rz(1.0452966) q[2];
sx q[2];
rz(-2.3901148) q[2];
sx q[2];
rz(0.95096248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54431984) q[1];
sx q[1];
rz(-1.2922799) q[1];
sx q[1];
rz(2.0803948) q[1];
rz(-0.3955807) q[3];
sx q[3];
rz(-2.1115344) q[3];
sx q[3];
rz(0.21025019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41844765) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(-2.1212228) q[2];
rz(2.4127035) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(-1.6093904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51434022) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(-0.038272055) q[0];
rz(0.12380869) q[1];
sx q[1];
rz(-0.47272155) q[1];
sx q[1];
rz(-1.5708057) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3338275) q[0];
sx q[0];
rz(-3.1196731) q[0];
sx q[0];
rz(-1.0967361) q[0];
rz(1.471631) q[2];
sx q[2];
rz(-0.80201521) q[2];
sx q[2];
rz(-0.85814171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96882183) q[1];
sx q[1];
rz(-1.5712196) q[1];
sx q[1];
rz(-1.569773) q[1];
rz(-pi) q[2];
rz(1.6917211) q[3];
sx q[3];
rz(-2.0837415) q[3];
sx q[3];
rz(2.4642022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6191285) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(0.5788571) q[2];
rz(1.045643) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(2.3845909) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4724562) q[0];
sx q[0];
rz(-1.0772912) q[0];
sx q[0];
rz(-0.45397595) q[0];
rz(0.16481915) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(1.4705315) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6413841) q[0];
sx q[0];
rz(-0.427503) q[0];
sx q[0];
rz(-1.1469202) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0963832) q[2];
sx q[2];
rz(-1.8586577) q[2];
sx q[2];
rz(-2.9381816) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3537703) q[1];
sx q[1];
rz(-0.90836891) q[1];
sx q[1];
rz(-1.1819928) q[1];
rz(-pi) q[2];
rz(0.65916797) q[3];
sx q[3];
rz(-1.1416832) q[3];
sx q[3];
rz(1.3339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7620324) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(1.3564159) q[2];
rz(-0.49946579) q[3];
sx q[3];
rz(-1.6127337) q[3];
sx q[3];
rz(-2.0904026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9143963) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(2.0830578) q[0];
rz(-1.9805485) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(-3.0243691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974221) q[0];
sx q[0];
rz(-1.6099842) q[0];
sx q[0];
rz(2.5990361) q[0];
rz(-pi) q[1];
rz(-2.9449844) q[2];
sx q[2];
rz(-2.2067382) q[2];
sx q[2];
rz(1.8064808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1438751) q[1];
sx q[1];
rz(-1.2923601) q[1];
sx q[1];
rz(-0.77634676) q[1];
x q[2];
rz(-2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(1.5269296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3566572) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(1.3680722) q[2];
rz(-1.28537) q[3];
sx q[3];
rz(-2.0338438) q[3];
sx q[3];
rz(-2.4135597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0990937) q[0];
sx q[0];
rz(-0.92785257) q[0];
sx q[0];
rz(1.7294783) q[0];
rz(1.0026576) q[1];
sx q[1];
rz(-0.24191562) q[1];
sx q[1];
rz(-1.5589421) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5413645) q[0];
sx q[0];
rz(-1.9892271) q[0];
sx q[0];
rz(-0.61516841) q[0];
x q[1];
rz(-1.6539668) q[2];
sx q[2];
rz(-2.2874333) q[2];
sx q[2];
rz(-2.4019474) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11410275) q[1];
sx q[1];
rz(-1.5943502) q[1];
sx q[1];
rz(2.1211758) q[1];
x q[2];
rz(3.047278) q[3];
sx q[3];
rz(-1.8911165) q[3];
sx q[3];
rz(-1.3195147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5530508) q[2];
sx q[2];
rz(-1.43247) q[2];
sx q[2];
rz(-2.7480965) q[2];
rz(0.95997512) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(0.967832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680962) q[0];
sx q[0];
rz(-3.0143026) q[0];
sx q[0];
rz(2.9582276) q[0];
rz(-0.49194899) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(1.2194182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26579612) q[0];
sx q[0];
rz(-1.6599433) q[0];
sx q[0];
rz(-1.6515948) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1111818) q[2];
sx q[2];
rz(-1.4923499) q[2];
sx q[2];
rz(0.65785656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1893434) q[1];
sx q[1];
rz(-0.5533411) q[1];
sx q[1];
rz(-2.4207741) q[1];
x q[2];
rz(1.7352571) q[3];
sx q[3];
rz(-2.6643848) q[3];
sx q[3];
rz(2.8754614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.4364028) q[2];
sx q[2];
rz(0.08237002) q[2];
rz(0.92665893) q[3];
sx q[3];
rz(-2.4121273) q[3];
sx q[3];
rz(-1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.6019186) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(-0.71640054) q[0];
rz(-0.40388233) q[1];
sx q[1];
rz(-1.5682181) q[1];
sx q[1];
rz(-1.7049888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2540909) q[0];
sx q[0];
rz(-1.9483101) q[0];
sx q[0];
rz(-2.98682) q[0];
rz(2.9173849) q[2];
sx q[2];
rz(-2.6022403) q[2];
sx q[2];
rz(1.2559079) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.389839) q[1];
sx q[1];
rz(-1.1734278) q[1];
sx q[1];
rz(-1.6786871) q[1];
rz(-pi) q[2];
rz(2.4190046) q[3];
sx q[3];
rz(-0.94959414) q[3];
sx q[3];
rz(-0.91739839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7776103) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(3.0762365) q[2];
rz(0.86723793) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6548178) q[0];
sx q[0];
rz(-1.4210533) q[0];
sx q[0];
rz(-2.4105893) q[0];
rz(-2.3664318) q[1];
sx q[1];
rz(-0.33439264) q[1];
sx q[1];
rz(2.7269272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054798445) q[0];
sx q[0];
rz(-0.86547422) q[0];
sx q[0];
rz(-2.2948625) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35279775) q[2];
sx q[2];
rz(-0.55856201) q[2];
sx q[2];
rz(-1.1272023) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.848661) q[1];
sx q[1];
rz(-2.7457016) q[1];
sx q[1];
rz(-0.69940523) q[1];
rz(-pi) q[2];
rz(2.25423) q[3];
sx q[3];
rz(-2.931224) q[3];
sx q[3];
rz(-0.31926814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23184648) q[2];
sx q[2];
rz(-0.48603386) q[2];
sx q[2];
rz(0.092378423) q[2];
rz(0.77629027) q[3];
sx q[3];
rz(-1.679136) q[3];
sx q[3];
rz(2.9379454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.48731503) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(0.017729433) q[0];
rz(1.0461294) q[1];
sx q[1];
rz(-1.7283864) q[1];
sx q[1];
rz(1.0606934) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720406) q[0];
sx q[0];
rz(-1.5020796) q[0];
sx q[0];
rz(1.5567354) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1938042) q[2];
sx q[2];
rz(-0.62219884) q[2];
sx q[2];
rz(0.54905984) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.880278) q[1];
sx q[1];
rz(-0.78431097) q[1];
sx q[1];
rz(-1.8159291) q[1];
rz(-0.10194998) q[3];
sx q[3];
rz(-0.52604874) q[3];
sx q[3];
rz(-0.40204429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1549418) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(-0.6443392) q[2];
rz(0.84195697) q[3];
sx q[3];
rz(-1.5717477) q[3];
sx q[3];
rz(-0.90698609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3802721) q[0];
sx q[0];
rz(-2.705882) q[0];
sx q[0];
rz(0.45183387) q[0];
rz(1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(-1.570328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1033882) q[0];
sx q[0];
rz(-1.4006091) q[0];
sx q[0];
rz(1.5645932) q[0];
rz(1.0958448) q[2];
sx q[2];
rz(-1.4670851) q[2];
sx q[2];
rz(-2.6748859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.328969) q[1];
sx q[1];
rz(-1.6996034) q[1];
sx q[1];
rz(1.179943) q[1];
x q[2];
rz(2.4191946) q[3];
sx q[3];
rz(-1.5649475) q[3];
sx q[3];
rz(1.7286144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62341225) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(1.5709467) q[2];
rz(3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-2.6251729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8138206) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(-2.9227921) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(2.0135638) q[2];
sx q[2];
rz(-1.762286) q[2];
sx q[2];
rz(1.0049835) q[2];
rz(-2.4841251) q[3];
sx q[3];
rz(-1.8651265) q[3];
sx q[3];
rz(-0.65262564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
