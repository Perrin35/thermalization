OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.83752051) q[0];
sx q[0];
rz(-0.66523319) q[0];
sx q[0];
rz(-1.9159777) q[0];
rz(-0.6839112) q[1];
sx q[1];
rz(3.5332503) q[1];
sx q[1];
rz(11.864301) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4944899) q[0];
sx q[0];
rz(-2.4328874) q[0];
sx q[0];
rz(-0.1811709) q[0];
rz(-1.4058787) q[2];
sx q[2];
rz(-1.5497297) q[2];
sx q[2];
rz(1.5677468) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4106284) q[1];
sx q[1];
rz(-2.0515039) q[1];
sx q[1];
rz(-0.60810535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8075712) q[3];
sx q[3];
rz(-0.70565685) q[3];
sx q[3];
rz(-1.2232085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2012653) q[2];
sx q[2];
rz(-1.5387115) q[2];
sx q[2];
rz(-2.8931457) q[2];
rz(-1.1343608) q[3];
sx q[3];
rz(-0.13392197) q[3];
sx q[3];
rz(1.1321446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0487173) q[0];
sx q[0];
rz(-0.55260125) q[0];
sx q[0];
rz(-1.5930814) q[0];
rz(0.50367194) q[1];
sx q[1];
rz(-1.9182938) q[1];
sx q[1];
rz(-2.7647387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6183313) q[0];
sx q[0];
rz(-2.0428162) q[0];
sx q[0];
rz(-1.5492155) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69825577) q[2];
sx q[2];
rz(-0.94329903) q[2];
sx q[2];
rz(2.1098532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.406817) q[1];
sx q[1];
rz(-1.8207014) q[1];
sx q[1];
rz(-0.49244778) q[1];
x q[2];
rz(0.21864076) q[3];
sx q[3];
rz(-2.8014516) q[3];
sx q[3];
rz(0.66375932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5571931) q[2];
sx q[2];
rz(-0.69584766) q[2];
sx q[2];
rz(2.2759571) q[2];
rz(-1.4731167) q[3];
sx q[3];
rz(-2.1560463) q[3];
sx q[3];
rz(2.1540811) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5675548) q[0];
sx q[0];
rz(-1.8542629) q[0];
sx q[0];
rz(-1.0512742) q[0];
rz(1.4569262) q[1];
sx q[1];
rz(-1.8345865) q[1];
sx q[1];
rz(1.6251224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3932289) q[0];
sx q[0];
rz(-0.81384515) q[0];
sx q[0];
rz(-0.39155829) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2411738) q[2];
sx q[2];
rz(-2.134765) q[2];
sx q[2];
rz(-1.054686) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.17371) q[1];
sx q[1];
rz(-0.62025242) q[1];
sx q[1];
rz(-0.8773251) q[1];
rz(-2.4184035) q[3];
sx q[3];
rz(-3.0602877) q[3];
sx q[3];
rz(-1.7618084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39525825) q[2];
sx q[2];
rz(-2.0863775) q[2];
sx q[2];
rz(2.9193817) q[2];
rz(2.639751) q[3];
sx q[3];
rz(-1.7343438) q[3];
sx q[3];
rz(-0.80880222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89268452) q[0];
sx q[0];
rz(-2.6539256) q[0];
sx q[0];
rz(-1.1647613) q[0];
rz(2.248863) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(-0.28183118) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.696071) q[0];
sx q[0];
rz(-1.0907249) q[0];
sx q[0];
rz(1.4351373) q[0];
rz(-2.149942) q[2];
sx q[2];
rz(-2.1985612) q[2];
sx q[2];
rz(2.6038458) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6960245) q[1];
sx q[1];
rz(-2.9640475) q[1];
sx q[1];
rz(-1.1648747) q[1];
x q[2];
rz(-0.025679703) q[3];
sx q[3];
rz(-1.1246944) q[3];
sx q[3];
rz(3.1334973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.07936) q[2];
sx q[2];
rz(-2.5793109) q[2];
sx q[2];
rz(-0.53323659) q[2];
rz(2.0594635) q[3];
sx q[3];
rz(-2.5366668) q[3];
sx q[3];
rz(-0.81348872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.767652) q[0];
sx q[0];
rz(-0.39893183) q[0];
sx q[0];
rz(-1.8178513) q[0];
rz(-0.81870493) q[1];
sx q[1];
rz(-2.3981514) q[1];
sx q[1];
rz(2.0735819) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6477953) q[0];
sx q[0];
rz(-1.4915221) q[0];
sx q[0];
rz(-2.981841) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3049561) q[2];
sx q[2];
rz(-1.0435487) q[2];
sx q[2];
rz(-3.1339147) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1966168) q[1];
sx q[1];
rz(-0.91755897) q[1];
sx q[1];
rz(1.7004299) q[1];
rz(-pi) q[2];
rz(2.0433918) q[3];
sx q[3];
rz(-1.7180746) q[3];
sx q[3];
rz(2.6781899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1800804) q[2];
sx q[2];
rz(-2.027812) q[2];
sx q[2];
rz(2.2799344) q[2];
rz(1.0263475) q[3];
sx q[3];
rz(-1.9828826) q[3];
sx q[3];
rz(1.1177184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.79477972) q[0];
sx q[0];
rz(-0.81951278) q[0];
sx q[0];
rz(-0.89282194) q[0];
rz(-0.18094856) q[1];
sx q[1];
rz(-2.4632958) q[1];
sx q[1];
rz(-0.13042626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8065211) q[0];
sx q[0];
rz(-1.4627883) q[0];
sx q[0];
rz(1.7301225) q[0];
rz(1.0662634) q[2];
sx q[2];
rz(-2.1158525) q[2];
sx q[2];
rz(-1.8503328) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6857026) q[1];
sx q[1];
rz(-1.8079385) q[1];
sx q[1];
rz(-0.030045998) q[1];
rz(-pi) q[2];
rz(0.1993066) q[3];
sx q[3];
rz(-0.79108566) q[3];
sx q[3];
rz(-1.4996604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9170407) q[2];
sx q[2];
rz(-1.9150534) q[2];
sx q[2];
rz(-2.2652333) q[2];
rz(-1.2419491) q[3];
sx q[3];
rz(-0.94049898) q[3];
sx q[3];
rz(2.131264) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4001615) q[0];
sx q[0];
rz(-0.09859666) q[0];
sx q[0];
rz(-2.7778991) q[0];
rz(0.74186507) q[1];
sx q[1];
rz(-1.0262841) q[1];
sx q[1];
rz(-2.738764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8217709) q[0];
sx q[0];
rz(-1.9812225) q[0];
sx q[0];
rz(-0.26217802) q[0];
rz(-pi) q[1];
rz(0.37781395) q[2];
sx q[2];
rz(-1.1973518) q[2];
sx q[2];
rz(2.6732973) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.070870474) q[1];
sx q[1];
rz(-1.6477343) q[1];
sx q[1];
rz(1.4489732) q[1];
rz(-1.2996583) q[3];
sx q[3];
rz(-0.95356546) q[3];
sx q[3];
rz(0.24880508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3591298) q[2];
sx q[2];
rz(-0.368258) q[2];
sx q[2];
rz(-1.5472319) q[2];
rz(2.3063229) q[3];
sx q[3];
rz(-0.95649496) q[3];
sx q[3];
rz(0.48867759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816724) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(0.51505995) q[0];
rz(1.7022279) q[1];
sx q[1];
rz(-1.2934338) q[1];
sx q[1];
rz(2.7424367) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8989152) q[0];
sx q[0];
rz(-1.5726046) q[0];
sx q[0];
rz(-0.13092069) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3005343) q[2];
sx q[2];
rz(-2.3944602) q[2];
sx q[2];
rz(-0.085937339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5134352) q[1];
sx q[1];
rz(-1.5062467) q[1];
sx q[1];
rz(-1.8316395) q[1];
rz(-pi) q[2];
rz(-0.065654556) q[3];
sx q[3];
rz(-1.5720417) q[3];
sx q[3];
rz(0.37855442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9814375) q[2];
sx q[2];
rz(-1.258054) q[2];
sx q[2];
rz(-1.0464279) q[2];
rz(2.4753172) q[3];
sx q[3];
rz(-0.25686887) q[3];
sx q[3];
rz(-2.090914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1995354) q[0];
sx q[0];
rz(-0.10534795) q[0];
sx q[0];
rz(-0.60607213) q[0];
rz(-2.3145158) q[1];
sx q[1];
rz(-1.8111633) q[1];
sx q[1];
rz(-1.8082632) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5510941) q[0];
sx q[0];
rz(-0.37305377) q[0];
sx q[0];
rz(-0.28753745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.016388) q[2];
sx q[2];
rz(-2.2632709) q[2];
sx q[2];
rz(-2.8852579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3279475) q[1];
sx q[1];
rz(-0.678181) q[1];
sx q[1];
rz(1.3840883) q[1];
rz(0.17154947) q[3];
sx q[3];
rz(-2.6772873) q[3];
sx q[3];
rz(-2.7363079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1313608) q[2];
sx q[2];
rz(-2.2183245) q[2];
sx q[2];
rz(-2.7566747) q[2];
rz(-0.63079232) q[3];
sx q[3];
rz(-1.2814949) q[3];
sx q[3];
rz(-1.7990641) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25403062) q[0];
sx q[0];
rz(-1.7429202) q[0];
sx q[0];
rz(1.4400462) q[0];
rz(-0.74132672) q[1];
sx q[1];
rz(-1.4011551) q[1];
sx q[1];
rz(2.8828566) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9567159) q[0];
sx q[0];
rz(-1.8167082) q[0];
sx q[0];
rz(-3.1374404) q[0];
rz(-2.0912295) q[2];
sx q[2];
rz(-0.58133091) q[2];
sx q[2];
rz(-2.5985825) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30686305) q[1];
sx q[1];
rz(-2.266699) q[1];
sx q[1];
rz(-0.45180288) q[1];
rz(-pi) q[2];
rz(-0.40269884) q[3];
sx q[3];
rz(-0.81953632) q[3];
sx q[3];
rz(1.0454901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4197293) q[2];
sx q[2];
rz(-2.0294919) q[2];
sx q[2];
rz(1.8224243) q[2];
rz(2.3223274) q[3];
sx q[3];
rz(-1.1300488) q[3];
sx q[3];
rz(2.5949902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7405613) q[0];
sx q[0];
rz(-0.87420976) q[0];
sx q[0];
rz(2.5921205) q[0];
rz(1.4389379) q[1];
sx q[1];
rz(-0.66822744) q[1];
sx q[1];
rz(-2.6034036) q[1];
rz(-0.86072154) q[2];
sx q[2];
rz(-1.5624983) q[2];
sx q[2];
rz(2.7539243) q[2];
rz(2.1323754) q[3];
sx q[3];
rz(-1.3604966) q[3];
sx q[3];
rz(0.58498912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
