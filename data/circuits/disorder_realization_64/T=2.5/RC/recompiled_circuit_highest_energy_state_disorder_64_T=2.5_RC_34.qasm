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
rz(-0.023802726) q[0];
sx q[0];
rz(4.2476141) q[0];
sx q[0];
rz(10.201693) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(-0.97631747) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7556954) q[0];
sx q[0];
rz(-1.160826) q[0];
sx q[0];
rz(-2.9958821) q[0];
x q[1];
rz(1.637865) q[2];
sx q[2];
rz(-1.0080976) q[2];
sx q[2];
rz(-0.85516155) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.088625535) q[1];
sx q[1];
rz(-1.8421116) q[1];
sx q[1];
rz(-2.0120828) q[1];
rz(-0.032194897) q[3];
sx q[3];
rz(-2.1436084) q[3];
sx q[3];
rz(-0.16530748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50600791) q[2];
sx q[2];
rz(-0.053746544) q[2];
sx q[2];
rz(-0.40679833) q[2];
rz(-2.9721416) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(-1.0725526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88103831) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(0.0037923092) q[0];
rz(-3.0637528) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(-2.8357764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7821976) q[0];
sx q[0];
rz(-0.54963028) q[0];
sx q[0];
rz(1.2122985) q[0];
rz(-pi) q[1];
rz(-0.75188101) q[2];
sx q[2];
rz(-1.6691748) q[2];
sx q[2];
rz(1.6461314) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9287323) q[1];
sx q[1];
rz(-1.1747735) q[1];
sx q[1];
rz(-2.3235882) q[1];
x q[2];
rz(1.5068077) q[3];
sx q[3];
rz(-1.631569) q[3];
sx q[3];
rz(1.2343386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1514312) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-2.4988417) q[2];
rz(-3.0691872) q[3];
sx q[3];
rz(-2.0824771) q[3];
sx q[3];
rz(-1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088242315) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(-0.15637583) q[0];
rz(0.0414255) q[1];
sx q[1];
rz(-0.62774575) q[1];
sx q[1];
rz(1.5511537) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0254733) q[0];
sx q[0];
rz(-1.0295233) q[0];
sx q[0];
rz(-1.3295637) q[0];
rz(-2.4855108) q[2];
sx q[2];
rz(-1.6245884) q[2];
sx q[2];
rz(2.057586) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.620365) q[1];
sx q[1];
rz(-1.6832192) q[1];
sx q[1];
rz(0.53037723) q[1];
rz(-0.82156397) q[3];
sx q[3];
rz(-0.85431495) q[3];
sx q[3];
rz(0.91148538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7330043) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(0.020922529) q[2];
rz(0.18713348) q[3];
sx q[3];
rz(-2.9360866) q[3];
sx q[3];
rz(-3.0278897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9631831) q[0];
sx q[0];
rz(-2.8091176) q[0];
sx q[0];
rz(-0.43854976) q[0];
rz(1.6167538) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(-2.8964892) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.995718) q[0];
sx q[0];
rz(-0.85332131) q[0];
sx q[0];
rz(0.11086734) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7510277) q[2];
sx q[2];
rz(-2.7465944) q[2];
sx q[2];
rz(0.67827889) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0673128) q[1];
sx q[1];
rz(-1.7099705) q[1];
sx q[1];
rz(2.1140631) q[1];
rz(-1.887759) q[3];
sx q[3];
rz(-0.84268236) q[3];
sx q[3];
rz(-1.1945981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(-2.8098246) q[2];
rz(2.6541384) q[3];
sx q[3];
rz(-1.0468227) q[3];
sx q[3];
rz(2.148518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29193923) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(-2.3680903) q[0];
rz(-1.9603112) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(1.7519417) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8256102) q[0];
sx q[0];
rz(-1.9709236) q[0];
sx q[0];
rz(1.9670301) q[0];
rz(-pi) q[1];
x q[1];
rz(2.543879) q[2];
sx q[2];
rz(-2.1897912) q[2];
sx q[2];
rz(-2.5368382) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2323245) q[1];
sx q[1];
rz(-1.0977912) q[1];
sx q[1];
rz(1.478341) q[1];
x q[2];
rz(-0.64719154) q[3];
sx q[3];
rz(-2.8464937) q[3];
sx q[3];
rz(-1.5403252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2191849) q[2];
sx q[2];
rz(-1.2097404) q[2];
sx q[2];
rz(-0.47214559) q[2];
rz(-1.2989429) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2209114) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(-2.8780908) q[0];
rz(2.0384516) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(0.37364328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81955302) q[0];
sx q[0];
rz(-0.094498903) q[0];
sx q[0];
rz(0.87847932) q[0];
rz(-pi) q[1];
rz(-2.2977826) q[2];
sx q[2];
rz(-1.6622953) q[2];
sx q[2];
rz(-1.3841656) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18446021) q[1];
sx q[1];
rz(-1.2533979) q[1];
sx q[1];
rz(-2.1613792) q[1];
x q[2];
rz(2.3832537) q[3];
sx q[3];
rz(-2.3985574) q[3];
sx q[3];
rz(0.61279994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0282447) q[2];
sx q[2];
rz(-0.17013203) q[2];
sx q[2];
rz(2.587758) q[2];
rz(1.3977741) q[3];
sx q[3];
rz(-0.60269409) q[3];
sx q[3];
rz(2.8288614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5572307) q[0];
sx q[0];
rz(-2.1350242) q[0];
sx q[0];
rz(-1.9118017) q[0];
rz(0.23904414) q[1];
sx q[1];
rz(-1.6300647) q[1];
sx q[1];
rz(-0.30034932) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2588033) q[0];
sx q[0];
rz(-2.9414294) q[0];
sx q[0];
rz(0.61997719) q[0];
rz(-pi) q[1];
rz(-1.0763758) q[2];
sx q[2];
rz(-0.9305312) q[2];
sx q[2];
rz(-0.51039052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21619851) q[1];
sx q[1];
rz(-1.5860737) q[1];
sx q[1];
rz(-0.00038860996) q[1];
rz(-pi) q[2];
rz(2.2836779) q[3];
sx q[3];
rz(-1.5893717) q[3];
sx q[3];
rz(-1.9481079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.131669) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(-0.24492502) q[2];
rz(-0.51982546) q[3];
sx q[3];
rz(-2.2782785) q[3];
sx q[3];
rz(0.68827099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7898665) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(-1.9785731) q[0];
rz(0.066896833) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(-1.012872) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6179919) q[0];
sx q[0];
rz(-0.80192425) q[0];
sx q[0];
rz(2.2973934) q[0];
rz(2.9628721) q[2];
sx q[2];
rz(-1.0105437) q[2];
sx q[2];
rz(2.5541039) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.098274577) q[1];
sx q[1];
rz(-1.2407082) q[1];
sx q[1];
rz(-1.9041474) q[1];
rz(-pi) q[2];
rz(2.6462534) q[3];
sx q[3];
rz(-0.74995774) q[3];
sx q[3];
rz(-1.0487674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11671994) q[2];
sx q[2];
rz(-0.97815424) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(-0.60574496) q[3];
sx q[3];
rz(-0.79234684) q[3];
sx q[3];
rz(2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.087273) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(-3.1291381) q[0];
rz(0.74673486) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(-0.27997231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6387647) q[0];
sx q[0];
rz(-1.6831213) q[0];
sx q[0];
rz(-3.0905484) q[0];
x q[1];
rz(-1.1367646) q[2];
sx q[2];
rz(-1.2723337) q[2];
sx q[2];
rz(3.0500183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38281967) q[1];
sx q[1];
rz(-2.3291322) q[1];
sx q[1];
rz(0.42940213) q[1];
rz(-pi) q[2];
rz(0.98484184) q[3];
sx q[3];
rz(-2.7142314) q[3];
sx q[3];
rz(-2.7387184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81194699) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(-2.9296056) q[2];
rz(-2.3181465) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14281808) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(-0.69277358) q[0];
rz(-0.57299262) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(-0.43100345) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8739024) q[0];
sx q[0];
rz(-0.55492102) q[0];
sx q[0];
rz(0.11124994) q[0];
rz(2.5415) q[2];
sx q[2];
rz(-1.1962435) q[2];
sx q[2];
rz(2.3698344) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1173874) q[1];
sx q[1];
rz(-2.3820602) q[1];
sx q[1];
rz(2.7717436) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0447356) q[3];
sx q[3];
rz(-2.5383679) q[3];
sx q[3];
rz(-0.67765498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7196322) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-2.5813622) q[2];
rz(2.6719921) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(2.4277021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2847168) q[0];
sx q[0];
rz(-1.4061883) q[0];
sx q[0];
rz(-1.0549369) q[0];
rz(0.84125413) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(-1.2011436) q[2];
sx q[2];
rz(-1.8664202) q[2];
sx q[2];
rz(2.0640434) q[2];
rz(-0.91184323) q[3];
sx q[3];
rz(-1.5494491) q[3];
sx q[3];
rz(-0.56964239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
