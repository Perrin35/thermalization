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
rz(0.12662521) q[0];
sx q[0];
rz(-1.5669444) q[0];
sx q[0];
rz(2.6259165) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(-0.74695865) q[1];
sx q[1];
rz(0.42626122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.212972) q[0];
sx q[0];
rz(-2.3317695) q[0];
sx q[0];
rz(-2.9024603) q[0];
x q[1];
rz(-1.8581763) q[2];
sx q[2];
rz(-2.0567908) q[2];
sx q[2];
rz(0.11746841) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30657739) q[1];
sx q[1];
rz(-1.6284429) q[1];
sx q[1];
rz(0.60204864) q[1];
rz(-2.8167546) q[3];
sx q[3];
rz(-1.0110613) q[3];
sx q[3];
rz(-1.4046275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8969741) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(2.8986325) q[2];
rz(-2.7729559) q[3];
sx q[3];
rz(-2.536085) q[3];
sx q[3];
rz(-1.9093556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31545562) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(-2.7742703) q[0];
rz(0.79633725) q[1];
sx q[1];
rz(-1.0581191) q[1];
sx q[1];
rz(-0.64250362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028051826) q[0];
sx q[0];
rz(-1.100825) q[0];
sx q[0];
rz(-2.0451106) q[0];
rz(-pi) q[1];
rz(-0.55134578) q[2];
sx q[2];
rz(-1.2228106) q[2];
sx q[2];
rz(2.6642193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31256235) q[1];
sx q[1];
rz(-1.0110657) q[1];
sx q[1];
rz(-2.1055431) q[1];
rz(-pi) q[2];
rz(-2.0525371) q[3];
sx q[3];
rz(-1.6828487) q[3];
sx q[3];
rz(-2.5555536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.502304) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(-1.7712234) q[2];
rz(0.227452) q[3];
sx q[3];
rz(-1.2496313) q[3];
sx q[3];
rz(-1.1915709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3962536) q[0];
sx q[0];
rz(-2.9662913) q[0];
sx q[0];
rz(-1.1981717) q[0];
rz(-1.0802957) q[1];
sx q[1];
rz(-0.21427576) q[1];
sx q[1];
rz(-3.0467196) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085233) q[0];
sx q[0];
rz(-2.108639) q[0];
sx q[0];
rz(0.30098029) q[0];
rz(-pi) q[1];
rz(1.0611141) q[2];
sx q[2];
rz(-0.72619263) q[2];
sx q[2];
rz(-0.92483556) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97040365) q[1];
sx q[1];
rz(-0.53817777) q[1];
sx q[1];
rz(1.6506877) q[1];
x q[2];
rz(-2.1839147) q[3];
sx q[3];
rz(-2.6827178) q[3];
sx q[3];
rz(-1.4754627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0723116) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(2.7351232) q[2];
rz(-1.1786002) q[3];
sx q[3];
rz(-1.7013763) q[3];
sx q[3];
rz(1.9788205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56226319) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(-0.33600268) q[0];
rz(2.621189) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(-0.10428183) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2662243) q[0];
sx q[0];
rz(-1.2422529) q[0];
sx q[0];
rz(-1.9383517) q[0];
x q[1];
rz(1.5849131) q[2];
sx q[2];
rz(-1.7272564) q[2];
sx q[2];
rz(-0.9597646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.1946009) q[1];
sx q[1];
rz(-0.90819383) q[1];
sx q[1];
rz(-2.6688982) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6882668) q[3];
sx q[3];
rz(-1.6764055) q[3];
sx q[3];
rz(-2.3345514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61863724) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(-1.9970419) q[2];
rz(-2.5942904) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(-1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.8892141) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(2.5872173) q[0];
rz(-2.2303708) q[1];
sx q[1];
rz(-1.8781885) q[1];
sx q[1];
rz(-2.8401781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9145935) q[0];
sx q[0];
rz(-2.1560139) q[0];
sx q[0];
rz(-0.27497681) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7784987) q[2];
sx q[2];
rz(-1.951521) q[2];
sx q[2];
rz(2.8247339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44228783) q[1];
sx q[1];
rz(-1.3584474) q[1];
sx q[1];
rz(-1.1816756) q[1];
x q[2];
rz(0.88121342) q[3];
sx q[3];
rz(-0.47166079) q[3];
sx q[3];
rz(0.72615964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1201799) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-3.1346698) q[2];
rz(0.93112469) q[3];
sx q[3];
rz(-1.1313063) q[3];
sx q[3];
rz(-1.1091703) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865006) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(-1.6078) q[0];
rz(-2.4389229) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(0.30219561) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55794278) q[0];
sx q[0];
rz(-1.4690225) q[0];
sx q[0];
rz(2.2324865) q[0];
rz(-pi) q[1];
rz(-0.40001656) q[2];
sx q[2];
rz(-0.75655327) q[2];
sx q[2];
rz(-3.0885901) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8540878) q[1];
sx q[1];
rz(-0.64784986) q[1];
sx q[1];
rz(-1.8604398) q[1];
rz(-0.052107776) q[3];
sx q[3];
rz(-1.4192389) q[3];
sx q[3];
rz(-1.0406756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2493784) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(0.92602473) q[2];
rz(-1.2146436) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(-0.27629575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72325426) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-2.0330698) q[0];
rz(-0.33379894) q[1];
sx q[1];
rz(-1.326694) q[1];
sx q[1];
rz(-2.0557859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.007581) q[0];
sx q[0];
rz(-1.7756878) q[0];
sx q[0];
rz(1.4896859) q[0];
x q[1];
rz(0.18098197) q[2];
sx q[2];
rz(-1.1753193) q[2];
sx q[2];
rz(2.4534695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9232417) q[1];
sx q[1];
rz(-2.3585547) q[1];
sx q[1];
rz(0.79773517) q[1];
rz(-pi) q[2];
rz(-2.2074039) q[3];
sx q[3];
rz(-1.459834) q[3];
sx q[3];
rz(-2.8257089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6256025) q[2];
sx q[2];
rz(-0.36449271) q[2];
sx q[2];
rz(1.9773352) q[2];
rz(2.8734251) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(2.3837762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6044354) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(1.1489768) q[0];
rz(-0.2941429) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(1.0424967) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9669078) q[0];
sx q[0];
rz(-2.3174441) q[0];
sx q[0];
rz(-2.0281726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22886054) q[2];
sx q[2];
rz(-0.80297744) q[2];
sx q[2];
rz(-0.91509089) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62952215) q[1];
sx q[1];
rz(-0.40813706) q[1];
sx q[1];
rz(0.44993181) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6076261) q[3];
sx q[3];
rz(-2.4772518) q[3];
sx q[3];
rz(-3.1337332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5440172) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(1.3647122) q[2];
rz(-0.65258604) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(-0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.515601) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(0.01734497) q[0];
rz(0.01677244) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(2.9877072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2543593) q[0];
sx q[0];
rz(-1.6463929) q[0];
sx q[0];
rz(-0.31799728) q[0];
rz(-pi) q[1];
rz(0.95005401) q[2];
sx q[2];
rz(-2.1851106) q[2];
sx q[2];
rz(-1.7194192) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9964136) q[1];
sx q[1];
rz(-2.6657678) q[1];
sx q[1];
rz(-2.4563172) q[1];
rz(-0.023970402) q[3];
sx q[3];
rz(-1.4949189) q[3];
sx q[3];
rz(-0.45437231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.953557) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(-2.6970862) q[2];
rz(2.9514173) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(-1.7688513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1403777) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(1.0706527) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.670198) q[1];
sx q[1];
rz(-0.79107034) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2915724) q[0];
sx q[0];
rz(-2.794461) q[0];
sx q[0];
rz(-1.3107383) q[0];
x q[1];
rz(-0.39761646) q[2];
sx q[2];
rz(-1.8194345) q[2];
sx q[2];
rz(1.719081) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2998878) q[1];
sx q[1];
rz(-1.1763402) q[1];
sx q[1];
rz(-1.9576555) q[1];
rz(-2.0759301) q[3];
sx q[3];
rz(-1.4544684) q[3];
sx q[3];
rz(1.1992421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4967686) q[2];
sx q[2];
rz(-0.18459979) q[2];
sx q[2];
rz(1.6688639) q[2];
rz(-1.5276927) q[3];
sx q[3];
rz(-2.0852641) q[3];
sx q[3];
rz(-1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2033757) q[0];
sx q[0];
rz(-1.4980409) q[0];
sx q[0];
rz(-0.71312755) q[0];
rz(-2.6279502) q[1];
sx q[1];
rz(-1.7084264) q[1];
sx q[1];
rz(1.4485566) q[1];
rz(-2.7265139) q[2];
sx q[2];
rz(-2.5998301) q[2];
sx q[2];
rz(2.0143581) q[2];
rz(1.1340101) q[3];
sx q[3];
rz(-0.81098771) q[3];
sx q[3];
rz(0.24843957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
