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
rz(5.5362267) q[1];
sx q[1];
rz(9.8510392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58881271) q[0];
sx q[0];
rz(-2.3511887) q[0];
sx q[0];
rz(-1.8145723) q[0];
x q[1];
rz(2.6381016) q[2];
sx q[2];
rz(-1.3174743) q[2];
sx q[2];
rz(-1.5905141) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3037881) q[1];
sx q[1];
rz(-0.96988867) q[1];
sx q[1];
rz(-1.5008885) q[1];
x q[2];
rz(2.1549752) q[3];
sx q[3];
rz(-1.8446577) q[3];
sx q[3];
rz(-0.010771839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2446186) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(-2.8986325) q[2];
rz(0.36863676) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(-1.2322371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31545562) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(2.7742703) q[0];
rz(-0.79633725) q[1];
sx q[1];
rz(-2.0834736) q[1];
sx q[1];
rz(2.499089) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2655661) q[0];
sx q[0];
rz(-0.6548223) q[0];
sx q[0];
rz(-0.73237082) q[0];
rz(-2.5359277) q[2];
sx q[2];
rz(-2.4993651) q[2];
sx q[2];
rz(-2.5544744) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31256235) q[1];
sx q[1];
rz(-1.0110657) q[1];
sx q[1];
rz(1.0360495) q[1];
x q[2];
rz(-1.3325464) q[3];
sx q[3];
rz(-2.6479911) q[3];
sx q[3];
rz(0.77405888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6392886) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(-1.7712234) q[2];
rz(-0.227452) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(-1.1915709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3962536) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(1.1981717) q[0];
rz(-2.0612969) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(0.094873039) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8137275) q[0];
sx q[0];
rz(-1.3133606) q[0];
sx q[0];
rz(1.0125005) q[0];
rz(0.91135613) q[2];
sx q[2];
rz(-1.9007287) q[2];
sx q[2];
rz(2.099769) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66902924) q[1];
sx q[1];
rz(-1.5298784) q[1];
sx q[1];
rz(-2.1075691) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8646144) q[3];
sx q[3];
rz(-1.9414895) q[3];
sx q[3];
rz(2.1408368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0723116) q[2];
sx q[2];
rz(-0.60439622) q[2];
sx q[2];
rz(0.4064694) q[2];
rz(-1.9629924) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.1627722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56226319) q[0];
sx q[0];
rz(-3.0070906) q[0];
sx q[0];
rz(0.33600268) q[0];
rz(-0.52040368) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(3.0373108) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18096237) q[0];
sx q[0];
rz(-1.2237566) q[0];
sx q[0];
rz(-0.35023679) q[0];
rz(2.9851172) q[2];
sx q[2];
rz(-1.5847407) q[2];
sx q[2];
rz(-2.5283611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49890624) q[1];
sx q[1];
rz(-2.3489174) q[1];
sx q[1];
rz(1.0426056) q[1];
rz(-0.4533259) q[3];
sx q[3];
rz(-1.4651872) q[3];
sx q[3];
rz(-2.3345514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61863724) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(1.9970419) q[2];
rz(-2.5942904) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(2.0539637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8892141) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(-0.55437535) q[0];
rz(2.2303708) q[1];
sx q[1];
rz(-1.8781885) q[1];
sx q[1];
rz(-0.30141452) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6431993) q[0];
sx q[0];
rz(-1.7991156) q[0];
sx q[0];
rz(-2.1737745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6658556) q[2];
sx q[2];
rz(-2.7103376) q[2];
sx q[2];
rz(-2.3088344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5382945) q[1];
sx q[1];
rz(-0.4406826) q[1];
sx q[1];
rz(2.0875817) q[1];
x q[2];
rz(0.88121342) q[3];
sx q[3];
rz(-2.6699319) q[3];
sx q[3];
rz(2.415433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.021412795) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-3.1346698) q[2];
rz(2.210468) q[3];
sx q[3];
rz(-1.1313063) q[3];
sx q[3];
rz(-2.0324223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329426) q[0];
sx q[0];
rz(-0.83704346) q[0];
sx q[0];
rz(-1.6078) q[0];
rz(-2.4389229) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(-2.839397) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9990197) q[0];
sx q[0];
rz(-0.66830615) q[0];
sx q[0];
rz(-1.4060941) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7415761) q[2];
sx q[2];
rz(-0.75655327) q[2];
sx q[2];
rz(3.0885901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.211765) q[1];
sx q[1];
rz(-2.1874839) q[1];
sx q[1];
rz(2.9287128) q[1];
rz(-pi) q[2];
rz(-1.8994529) q[3];
sx q[3];
rz(-2.9813926) q[3];
sx q[3];
rz(2.4335308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.89221421) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(2.2155679) q[2];
rz(1.9269491) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(2.8652969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4183384) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(-2.8077937) q[1];
sx q[1];
rz(-1.326694) q[1];
sx q[1];
rz(2.0557859) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1340116) q[0];
sx q[0];
rz(-1.3659048) q[0];
sx q[0];
rz(-1.4896859) q[0];
rz(-pi) q[1];
x q[1];
rz(1.169431) q[2];
sx q[2];
rz(-1.4039206) q[2];
sx q[2];
rz(-2.3292975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2183509) q[1];
sx q[1];
rz(-2.3585547) q[1];
sx q[1];
rz(-0.79773517) q[1];
rz(-pi) q[2];
rz(-2.2074039) q[3];
sx q[3];
rz(-1.6817586) q[3];
sx q[3];
rz(-0.31588376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51599017) q[2];
sx q[2];
rz(-2.7770999) q[2];
sx q[2];
rz(-1.1642574) q[2];
rz(-2.8734251) q[3];
sx q[3];
rz(-1.8987013) q[3];
sx q[3];
rz(2.3837762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5371573) q[0];
sx q[0];
rz(-0.51979655) q[0];
sx q[0];
rz(-1.9926158) q[0];
rz(-0.2941429) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(1.0424967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54759083) q[0];
sx q[0];
rz(-2.2896575) q[0];
sx q[0];
rz(-2.6963364) q[0];
rz(-pi) q[1];
rz(2.3518219) q[2];
sx q[2];
rz(-1.7347448) q[2];
sx q[2];
rz(2.3254834) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62952215) q[1];
sx q[1];
rz(-0.40813706) q[1];
sx q[1];
rz(0.44993181) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90678471) q[3];
sx q[3];
rz(-1.5480925) q[3];
sx q[3];
rz(-1.5496538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59757549) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(-1.7768804) q[2];
rz(0.65258604) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(-2.4932388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.515601) q[0];
sx q[0];
rz(-2.7643272) q[0];
sx q[0];
rz(0.01734497) q[0];
rz(0.01677244) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(2.9877072) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8003004) q[0];
sx q[0];
rz(-1.8878536) q[0];
sx q[0];
rz(1.6503667) q[0];
rz(-2.4520313) q[2];
sx q[2];
rz(-0.84367263) q[2];
sx q[2];
rz(-0.52983701) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20275234) q[1];
sx q[1];
rz(-1.2766663) q[1];
sx q[1];
rz(0.37962706) q[1];
x q[2];
rz(3.1176223) q[3];
sx q[3];
rz(-1.4949189) q[3];
sx q[3];
rz(2.6872203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.953557) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(-2.6970862) q[2];
rz(2.9514173) q[3];
sx q[3];
rz(-1.5835652) q[3];
sx q[3];
rz(-1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.001215) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(-1.0706527) q[0];
rz(-0.045914687) q[1];
sx q[1];
rz(-1.670198) q[1];
sx q[1];
rz(2.3505223) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.175638) q[0];
sx q[0];
rz(-1.6583867) q[0];
sx q[0];
rz(1.2344633) q[0];
rz(-pi) q[1];
rz(-2.7439762) q[2];
sx q[2];
rz(-1.3221581) q[2];
sx q[2];
rz(1.719081) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4262096) q[1];
sx q[1];
rz(-1.2150303) q[1];
sx q[1];
rz(-0.4224311) q[1];
rz(-2.0759301) q[3];
sx q[3];
rz(-1.6871243) q[3];
sx q[3];
rz(-1.1992421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4967686) q[2];
sx q[2];
rz(-2.9569929) q[2];
sx q[2];
rz(-1.4727288) q[2];
rz(1.5276927) q[3];
sx q[3];
rz(-1.0563285) q[3];
sx q[3];
rz(1.9014026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.938217) q[0];
sx q[0];
rz(-1.6435517) q[0];
sx q[0];
rz(2.4284651) q[0];
rz(-2.6279502) q[1];
sx q[1];
rz(-1.7084264) q[1];
sx q[1];
rz(1.4485566) q[1];
rz(1.8088874) q[2];
sx q[2];
rz(-1.0793964) q[2];
sx q[2];
rz(2.4894077) q[2];
rz(-0.80908262) q[3];
sx q[3];
rz(-1.8825023) q[3];
sx q[3];
rz(-1.0112345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
