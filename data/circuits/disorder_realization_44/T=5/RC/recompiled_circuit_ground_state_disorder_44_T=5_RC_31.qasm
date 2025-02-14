OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(-1.4491117) q[0];
sx q[0];
rz(1.4504855) q[0];
rz(-2.1679572) q[1];
sx q[1];
rz(-1.4373625) q[1];
sx q[1];
rz(0.91926423) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5333119) q[0];
sx q[0];
rz(-1.5402147) q[0];
sx q[0];
rz(1.5879059) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4798574) q[2];
sx q[2];
rz(-1.6906066) q[2];
sx q[2];
rz(2.7594942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3340993) q[1];
sx q[1];
rz(-1.4791282) q[1];
sx q[1];
rz(-0.86218545) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1000865) q[3];
sx q[3];
rz(-1.19095) q[3];
sx q[3];
rz(-2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8093402) q[2];
sx q[2];
rz(-1.4490178) q[2];
sx q[2];
rz(1.7285041) q[2];
rz(0.20279065) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(-0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.24421144) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(0.71075034) q[0];
rz(2.3731025) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(-1.01952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43895129) q[0];
sx q[0];
rz(-2.0176689) q[0];
sx q[0];
rz(1.5269482) q[0];
x q[1];
rz(-1.7350082) q[2];
sx q[2];
rz(-0.43638849) q[2];
sx q[2];
rz(-1.1498888) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72538917) q[1];
sx q[1];
rz(-0.37230154) q[1];
sx q[1];
rz(-1.6561693) q[1];
rz(-pi) q[2];
rz(-2.1250167) q[3];
sx q[3];
rz(-1.4999522) q[3];
sx q[3];
rz(1.4431825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3780313) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(1.2126806) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0843622) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(-2.6257264) q[1];
sx q[1];
rz(-0.54793826) q[1];
sx q[1];
rz(-0.9224433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8144174) q[0];
sx q[0];
rz(-1.6285768) q[0];
sx q[0];
rz(-1.5301276) q[0];
rz(2.0751245) q[2];
sx q[2];
rz(-2.2309003) q[2];
sx q[2];
rz(-0.65048993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53239142) q[1];
sx q[1];
rz(-0.89840404) q[1];
sx q[1];
rz(-2.3228541) q[1];
x q[2];
rz(-2.8124468) q[3];
sx q[3];
rz(-1.7308047) q[3];
sx q[3];
rz(3.126006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.018365232) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(1.3398735) q[2];
rz(-2.5943622) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(-2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0860586) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(2.9823629) q[0];
rz(0.010146443) q[1];
sx q[1];
rz(-1.0228913) q[1];
sx q[1];
rz(2.8841282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2672294) q[0];
sx q[0];
rz(-0.82634514) q[0];
sx q[0];
rz(-2.6494725) q[0];
rz(-3.0765216) q[2];
sx q[2];
rz(-2.3461968) q[2];
sx q[2];
rz(0.8000904) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83492571) q[1];
sx q[1];
rz(-1.67027) q[1];
sx q[1];
rz(-2.5132781) q[1];
x q[2];
rz(-0.80735029) q[3];
sx q[3];
rz(-2.0424543) q[3];
sx q[3];
rz(1.1273257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(-2.6687458) q[2];
rz(2.4371448) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5064297) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(0.065486431) q[0];
rz(-0.40924117) q[1];
sx q[1];
rz(-1.1301273) q[1];
sx q[1];
rz(-1.7154891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62140853) q[0];
sx q[0];
rz(-1.6637633) q[0];
sx q[0];
rz(0.21958242) q[0];
x q[1];
rz(-0.3887811) q[2];
sx q[2];
rz(-1.6531303) q[2];
sx q[2];
rz(1.3363234) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8219442) q[1];
sx q[1];
rz(-1.2927755) q[1];
sx q[1];
rz(-1.3389498) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0423921) q[3];
sx q[3];
rz(-0.33337731) q[3];
sx q[3];
rz(2.7957145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19491974) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(1.8348414) q[2];
rz(-1.9715747) q[3];
sx q[3];
rz(-0.40478671) q[3];
sx q[3];
rz(-3.034333) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65866798) q[0];
sx q[0];
rz(-1.985745) q[0];
sx q[0];
rz(-2.7401127) q[0];
rz(2.4688156) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(-2.2875517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14388785) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(-1.7884939) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18547345) q[2];
sx q[2];
rz(-1.4234241) q[2];
sx q[2];
rz(-2.3285248) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6352977) q[1];
sx q[1];
rz(-2.1599401) q[1];
sx q[1];
rz(1.4511257) q[1];
x q[2];
rz(0.60131945) q[3];
sx q[3];
rz(-2.4375705) q[3];
sx q[3];
rz(2.6773334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65471571) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(2.5901637) q[3];
sx q[3];
rz(-1.5588372) q[3];
sx q[3];
rz(0.83561713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42171445) q[0];
sx q[0];
rz(-0.28110176) q[0];
sx q[0];
rz(1.9792492) q[0];
rz(1.989919) q[1];
sx q[1];
rz(-1.78777) q[1];
sx q[1];
rz(0.91845671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29354039) q[0];
sx q[0];
rz(-1.2402724) q[0];
sx q[0];
rz(2.8577515) q[0];
x q[1];
rz(-1.880391) q[2];
sx q[2];
rz(-2.4221276) q[2];
sx q[2];
rz(-1.3165064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20275252) q[1];
sx q[1];
rz(-1.5781286) q[1];
sx q[1];
rz(-1.5775024) q[1];
rz(-pi) q[2];
rz(0.17042589) q[3];
sx q[3];
rz(-1.9099351) q[3];
sx q[3];
rz(-0.91811524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.034417001) q[2];
sx q[2];
rz(-1.3849266) q[2];
sx q[2];
rz(2.6386063) q[2];
rz(0.77477396) q[3];
sx q[3];
rz(-2.6796902) q[3];
sx q[3];
rz(-1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4485432) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.6402798) q[0];
rz(0.34128183) q[1];
sx q[1];
rz(-1.7106067) q[1];
sx q[1];
rz(-1.1192082) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687896) q[0];
sx q[0];
rz(-0.69376341) q[0];
sx q[0];
rz(2.540178) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0430492) q[2];
sx q[2];
rz(-2.611428) q[2];
sx q[2];
rz(1.4250172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95204845) q[1];
sx q[1];
rz(-0.91390007) q[1];
sx q[1];
rz(1.7236962) q[1];
rz(-pi) q[2];
rz(-2.8265727) q[3];
sx q[3];
rz(-2.7286988) q[3];
sx q[3];
rz(1.1668432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1559653) q[2];
sx q[2];
rz(-0.95981821) q[2];
sx q[2];
rz(0.36925527) q[2];
rz(-0.30820942) q[3];
sx q[3];
rz(-2.1741368) q[3];
sx q[3];
rz(-1.6109899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592598) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(-2.7985213) q[0];
rz(1.169091) q[1];
sx q[1];
rz(-1.3645376) q[1];
sx q[1];
rz(-1.4287359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0562623) q[0];
sx q[0];
rz(-1.3199727) q[0];
sx q[0];
rz(2.73231) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4506017) q[2];
sx q[2];
rz(-0.96110247) q[2];
sx q[2];
rz(0.26870773) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45537585) q[1];
sx q[1];
rz(-1.1955402) q[1];
sx q[1];
rz(-0.065836716) q[1];
rz(-pi) q[2];
rz(-1.021528) q[3];
sx q[3];
rz(-0.77199751) q[3];
sx q[3];
rz(-0.87024161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(1.3915871) q[2];
rz(0.40677795) q[3];
sx q[3];
rz(-1.4429561) q[3];
sx q[3];
rz(2.2627635) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0402891) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(1.3579177) q[0];
rz(-1.9039924) q[1];
sx q[1];
rz(-1.0126746) q[1];
sx q[1];
rz(2.3416669) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3615204) q[0];
sx q[0];
rz(-1.8006386) q[0];
sx q[0];
rz(-0.80432463) q[0];
x q[1];
rz(1.2679891) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(-2.985266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.25617304) q[1];
sx q[1];
rz(-2.2347663) q[1];
sx q[1];
rz(1.7682942) q[1];
rz(2.8915845) q[3];
sx q[3];
rz(-2.0431314) q[3];
sx q[3];
rz(-1.3046622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7633729) q[2];
sx q[2];
rz(-1.4241445) q[2];
sx q[2];
rz(3.0685032) q[2];
rz(0.25589219) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(1.6528486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.273461) q[0];
sx q[0];
rz(-1.4556226) q[0];
sx q[0];
rz(1.8709394) q[0];
rz(-1.801626) q[1];
sx q[1];
rz(-1.3506964) q[1];
sx q[1];
rz(0.66257308) q[1];
rz(0.70941464) q[2];
sx q[2];
rz(-1.5723036) q[2];
sx q[2];
rz(0.2607762) q[2];
rz(-0.70684915) q[3];
sx q[3];
rz(-0.77597386) q[3];
sx q[3];
rz(1.5472277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
