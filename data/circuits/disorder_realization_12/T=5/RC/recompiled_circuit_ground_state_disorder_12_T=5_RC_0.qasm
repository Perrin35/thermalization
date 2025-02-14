OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(-1.8422814) q[0];
sx q[0];
rz(0.045825034) q[0];
rz(1.2996281) q[1];
sx q[1];
rz(-1.8027432) q[1];
sx q[1];
rz(-1.8884698) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6413413) q[0];
sx q[0];
rz(-1.1317028) q[0];
sx q[0];
rz(0.27089775) q[0];
rz(-pi) q[1];
rz(1.6145176) q[2];
sx q[2];
rz(-1.3695903) q[2];
sx q[2];
rz(-1.218856) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3435077) q[1];
sx q[1];
rz(-1.6899334) q[1];
sx q[1];
rz(0.75372523) q[1];
rz(-1.6016989) q[3];
sx q[3];
rz(-0.12750827) q[3];
sx q[3];
rz(0.37653437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32156285) q[2];
sx q[2];
rz(-1.4048615) q[2];
sx q[2];
rz(-0.31537867) q[2];
rz(0.086221181) q[3];
sx q[3];
rz(-3.0486221) q[3];
sx q[3];
rz(0.84403795) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.740199) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(0.33541086) q[0];
rz(-2.7566578) q[1];
sx q[1];
rz(-2.3454869) q[1];
sx q[1];
rz(-1.8523432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1803363) q[0];
sx q[0];
rz(-1.5085876) q[0];
sx q[0];
rz(-0.38952413) q[0];
rz(-pi) q[1];
rz(1.697549) q[2];
sx q[2];
rz(-0.64536649) q[2];
sx q[2];
rz(-0.54189787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3284387) q[1];
sx q[1];
rz(-2.6578609) q[1];
sx q[1];
rz(2.4991922) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3277169) q[3];
sx q[3];
rz(-1.4395003) q[3];
sx q[3];
rz(0.32570244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59911597) q[2];
sx q[2];
rz(-1.7581538) q[2];
sx q[2];
rz(1.4625589) q[2];
rz(-2.2595432) q[3];
sx q[3];
rz(-2.4817011) q[3];
sx q[3];
rz(-1.6409469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7733234) q[0];
sx q[0];
rz(-2.6550846) q[0];
sx q[0];
rz(-0.64273709) q[0];
rz(0.87359387) q[1];
sx q[1];
rz(-1.3106376) q[1];
sx q[1];
rz(-0.98683039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9844279) q[0];
sx q[0];
rz(-1.4207977) q[0];
sx q[0];
rz(-0.87777975) q[0];
rz(-pi) q[1];
rz(3.0796649) q[2];
sx q[2];
rz(-1.0871097) q[2];
sx q[2];
rz(2.7360423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1912811) q[1];
sx q[1];
rz(-2.1149968) q[1];
sx q[1];
rz(-2.4978375) q[1];
x q[2];
rz(-0.049792265) q[3];
sx q[3];
rz(-1.9450499) q[3];
sx q[3];
rz(0.58993898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70283908) q[2];
sx q[2];
rz(-3.1162016) q[2];
sx q[2];
rz(1.0566443) q[2];
rz(-1.3422525) q[3];
sx q[3];
rz(-1.4547576) q[3];
sx q[3];
rz(-0.87856436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.037828) q[0];
sx q[0];
rz(-0.28452888) q[0];
sx q[0];
rz(-1.1430662) q[0];
rz(0.68483886) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(-0.24982223) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9651325) q[0];
sx q[0];
rz(-1.1145381) q[0];
sx q[0];
rz(2.3224823) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5742214) q[2];
sx q[2];
rz(-1.5674233) q[2];
sx q[2];
rz(0.33004883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64678363) q[1];
sx q[1];
rz(-2.3338291) q[1];
sx q[1];
rz(-0.69112063) q[1];
rz(-pi) q[2];
rz(-0.0015578702) q[3];
sx q[3];
rz(-1.8127725) q[3];
sx q[3];
rz(0.93780692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41986856) q[2];
sx q[2];
rz(-1.3519752) q[2];
sx q[2];
rz(1.295759) q[2];
rz(1.3755679) q[3];
sx q[3];
rz(-1.7121168) q[3];
sx q[3];
rz(0.099055722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38882035) q[0];
sx q[0];
rz(-1.7701912) q[0];
sx q[0];
rz(-0.8465299) q[0];
rz(1.9580152) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(1.902098) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13279937) q[0];
sx q[0];
rz(-1.4504878) q[0];
sx q[0];
rz(-1.3161206) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6651004) q[2];
sx q[2];
rz(-0.38151729) q[2];
sx q[2];
rz(-0.028057773) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.073055) q[1];
sx q[1];
rz(-3.0843711) q[1];
sx q[1];
rz(-1.3609318) q[1];
rz(2.1679513) q[3];
sx q[3];
rz(-1.5212562) q[3];
sx q[3];
rz(1.3508391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5014629) q[2];
sx q[2];
rz(-1.5093466) q[2];
sx q[2];
rz(0.35422361) q[2];
rz(0.7192449) q[3];
sx q[3];
rz(-2.5651599) q[3];
sx q[3];
rz(0.80938068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(1.7648014) q[0];
sx q[0];
rz(-2.9053595) q[0];
sx q[0];
rz(0.34673044) q[0];
rz(-0.72225839) q[1];
sx q[1];
rz(-1.5303231) q[1];
sx q[1];
rz(-0.17682704) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8943902) q[0];
sx q[0];
rz(-0.48485928) q[0];
sx q[0];
rz(1.8762183) q[0];
rz(1.0718623) q[2];
sx q[2];
rz(-1.4889354) q[2];
sx q[2];
rz(0.53033842) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8521148) q[1];
sx q[1];
rz(-1.3657709) q[1];
sx q[1];
rz(-1.6774898) q[1];
rz(-2.1336055) q[3];
sx q[3];
rz(-1.6567818) q[3];
sx q[3];
rz(0.028349625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.56883812) q[2];
sx q[2];
rz(-0.73358959) q[2];
sx q[2];
rz(1.0432165) q[2];
rz(-2.9853232) q[3];
sx q[3];
rz(-2.0782491) q[3];
sx q[3];
rz(0.094303057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909376) q[0];
sx q[0];
rz(-2.1112878) q[0];
sx q[0];
rz(-1.913273) q[0];
rz(0.43371513) q[1];
sx q[1];
rz(-0.80077306) q[1];
sx q[1];
rz(2.0965651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997847) q[0];
sx q[0];
rz(-1.5130454) q[0];
sx q[0];
rz(-1.367913) q[0];
rz(-pi) q[1];
rz(-2.074998) q[2];
sx q[2];
rz(-2.3269751) q[2];
sx q[2];
rz(-2.2058918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4853128) q[1];
sx q[1];
rz(-1.6846141) q[1];
sx q[1];
rz(0.95358221) q[1];
x q[2];
rz(1.0492453) q[3];
sx q[3];
rz(-0.86498596) q[3];
sx q[3];
rz(-1.7569831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9016483) q[2];
sx q[2];
rz(-1.0239235) q[2];
sx q[2];
rz(-0.18590064) q[2];
rz(1.2402041) q[3];
sx q[3];
rz(-1.5408206) q[3];
sx q[3];
rz(-2.127229) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1321201) q[0];
sx q[0];
rz(-1.8931696) q[0];
sx q[0];
rz(0.15952071) q[0];
rz(-2.3347143) q[1];
sx q[1];
rz(-1.9069549) q[1];
sx q[1];
rz(-3.1330915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2293866) q[0];
sx q[0];
rz(-1.9766146) q[0];
sx q[0];
rz(-2.9738369) q[0];
rz(-pi) q[1];
rz(0.10430704) q[2];
sx q[2];
rz(-2.3416714) q[2];
sx q[2];
rz(1.6751573) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3169049) q[1];
sx q[1];
rz(-1.3456235) q[1];
sx q[1];
rz(-1.4359739) q[1];
rz(2.9656677) q[3];
sx q[3];
rz(-2.0485544) q[3];
sx q[3];
rz(-2.0328409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8622417) q[2];
sx q[2];
rz(-1.778435) q[2];
sx q[2];
rz(2.9478493) q[2];
rz(-1.6019609) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(1.1295454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2368471) q[0];
sx q[0];
rz(-1.3847677) q[0];
sx q[0];
rz(-0.73614502) q[0];
rz(-0.88788095) q[1];
sx q[1];
rz(-0.45279756) q[1];
sx q[1];
rz(0.050051659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0808635) q[0];
sx q[0];
rz(-1.2300778) q[0];
sx q[0];
rz(-1.3843892) q[0];
rz(-pi) q[1];
rz(2.8269269) q[2];
sx q[2];
rz(-1.6648714) q[2];
sx q[2];
rz(1.5794181) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60683317) q[1];
sx q[1];
rz(-1.8799361) q[1];
sx q[1];
rz(1.7844132) q[1];
rz(-pi) q[2];
rz(2.3042559) q[3];
sx q[3];
rz(-0.76238197) q[3];
sx q[3];
rz(1.1632944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8686409) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(0.0090553332) q[2];
rz(-2.791259) q[3];
sx q[3];
rz(-2.83559) q[3];
sx q[3];
rz(2.8370324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7250799) q[0];
sx q[0];
rz(-2.0820936) q[0];
sx q[0];
rz(0.9730202) q[0];
rz(-1.5618207) q[1];
sx q[1];
rz(-0.90619722) q[1];
sx q[1];
rz(2.3320893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0334629) q[0];
sx q[0];
rz(-1.9282883) q[0];
sx q[0];
rz(-0.78127677) q[0];
rz(-pi) q[1];
rz(-1.3415178) q[2];
sx q[2];
rz(-1.3827168) q[2];
sx q[2];
rz(3.0517202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36724801) q[1];
sx q[1];
rz(-1.4092236) q[1];
sx q[1];
rz(-2.0366686) q[1];
rz(2.1511492) q[3];
sx q[3];
rz(-1.1709612) q[3];
sx q[3];
rz(2.4357908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1805264) q[2];
sx q[2];
rz(-2.3193391) q[2];
sx q[2];
rz(-2.6194438) q[2];
rz(0.3565878) q[3];
sx q[3];
rz(-2.5643189) q[3];
sx q[3];
rz(3.0780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6575573) q[0];
sx q[0];
rz(-0.35150305) q[0];
sx q[0];
rz(1.2765314) q[0];
rz(0.36318489) q[1];
sx q[1];
rz(-0.52774944) q[1];
sx q[1];
rz(-1.48762) q[1];
rz(1.034201) q[2];
sx q[2];
rz(-2.2727439) q[2];
sx q[2];
rz(-1.5287413) q[2];
rz(-0.57989953) q[3];
sx q[3];
rz(-1.2968466) q[3];
sx q[3];
rz(-2.235835) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
