OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9085812) q[0];
sx q[0];
rz(-1.9549978) q[0];
sx q[0];
rz(-1.1458122) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(-0.29247984) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777893) q[0];
sx q[0];
rz(-1.6262494) q[0];
sx q[0];
rz(1.9652912) q[0];
rz(-pi) q[1];
rz(-2.3621759) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(1.2365637) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7293207) q[1];
sx q[1];
rz(-0.97727697) q[1];
sx q[1];
rz(-0.055298474) q[1];
x q[2];
rz(1.5745893) q[3];
sx q[3];
rz(-0.47382254) q[3];
sx q[3];
rz(-2.4592318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4108489) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(-0.4494108) q[2];
rz(-2.4959026) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(-1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21355024) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(-0.74203062) q[0];
rz(-1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.3630294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.610696) q[0];
sx q[0];
rz(-0.96141059) q[0];
sx q[0];
rz(1.0884398) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7861373) q[2];
sx q[2];
rz(-2.0305579) q[2];
sx q[2];
rz(-2.721867) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7098397) q[1];
sx q[1];
rz(-1.5515986) q[1];
sx q[1];
rz(-2.8019816) q[1];
x q[2];
rz(0.19248776) q[3];
sx q[3];
rz(-1.3093595) q[3];
sx q[3];
rz(0.9769494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(1.2724686) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(-3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0740046) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(0.2628251) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(1.2353108) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10997406) q[0];
sx q[0];
rz(-1.5556766) q[0];
sx q[0];
rz(-1.5605687) q[0];
x q[1];
rz(-0.026651816) q[2];
sx q[2];
rz(-1.3304552) q[2];
sx q[2];
rz(-1.196256) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3706627) q[1];
sx q[1];
rz(-1.2491033) q[1];
sx q[1];
rz(-1.6817723) q[1];
rz(-pi) q[2];
rz(1.9800817) q[3];
sx q[3];
rz(-2.403879) q[3];
sx q[3];
rz(2.0091332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(2.9612605) q[2];
rz(2.5056433) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(0.71907991) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(-1.0346574) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9229139) q[0];
sx q[0];
rz(-1.551034) q[0];
sx q[0];
rz(-3.0826871) q[0];
rz(1.7059533) q[2];
sx q[2];
rz(-1.5237336) q[2];
sx q[2];
rz(-1.0263718) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3116341) q[1];
sx q[1];
rz(-1.2882075) q[1];
sx q[1];
rz(2.4294873) q[1];
rz(0.33200522) q[3];
sx q[3];
rz(-2.6453291) q[3];
sx q[3];
rz(-0.23976025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93306142) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(0.92102712) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-2.5193118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3233933) q[0];
sx q[0];
rz(-0.33938956) q[0];
sx q[0];
rz(-2.1470977) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4786517) q[2];
sx q[2];
rz(-2.2170057) q[2];
sx q[2];
rz(-2.9421633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29187782) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(0.090976322) q[1];
rz(-pi) q[2];
rz(0.673224) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(2.4622963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.92419147) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(0.021281555) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(-2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(-3.0969627) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-3.1071641) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8196626) q[0];
sx q[0];
rz(-0.18580431) q[0];
sx q[0];
rz(1.6447322) q[0];
rz(-1.351864) q[2];
sx q[2];
rz(-2.6282675) q[2];
sx q[2];
rz(2.9190612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.039636314) q[1];
sx q[1];
rz(-0.64412457) q[1];
sx q[1];
rz(-2.1307751) q[1];
x q[2];
rz(1.4061635) q[3];
sx q[3];
rz(-0.87222404) q[3];
sx q[3];
rz(-2.5005831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3520711) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(2.1748523) q[2];
rz(-3.1294075) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(-0.99408856) q[0];
rz(-3.0864691) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(-0.73928839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2544884) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(-1.4332921) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67989345) q[2];
sx q[2];
rz(-0.69603622) q[2];
sx q[2];
rz(1.807715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97923219) q[1];
sx q[1];
rz(-1.7602966) q[1];
sx q[1];
rz(-0.53342553) q[1];
rz(0.61974157) q[3];
sx q[3];
rz(-1.3337787) q[3];
sx q[3];
rz(-1.851351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(2.901315) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(-1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(2.5575496) q[0];
rz(-2.7208327) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(-0.27841321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4415092) q[0];
sx q[0];
rz(-2.1372876) q[0];
sx q[0];
rz(-0.14871116) q[0];
rz(-pi) q[1];
rz(2.7483447) q[2];
sx q[2];
rz(-1.7448145) q[2];
sx q[2];
rz(0.72088036) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1380996) q[1];
sx q[1];
rz(-1.7894735) q[1];
sx q[1];
rz(-1.2117282) q[1];
rz(-1.9592459) q[3];
sx q[3];
rz(-2.5591345) q[3];
sx q[3];
rz(0.1633446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86928308) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(-2.753624) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(-0.2750245) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-2.9751119) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(-3.0678715) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8548944) q[0];
sx q[0];
rz(-0.44730967) q[0];
sx q[0];
rz(1.3351424) q[0];
rz(-0.080967112) q[2];
sx q[2];
rz(-1.4514187) q[2];
sx q[2];
rz(0.93842426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8718865) q[1];
sx q[1];
rz(-1.5333813) q[1];
sx q[1];
rz(2.5579631) q[1];
rz(-1.4599667) q[3];
sx q[3];
rz(-1.0004527) q[3];
sx q[3];
rz(3.0487206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5471197) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(0.42868844) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(1.930442) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4093032) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(1.0428585) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(1.7932549) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52085984) q[0];
sx q[0];
rz(-0.11611406) q[0];
sx q[0];
rz(1.7945047) q[0];
rz(-pi) q[1];
rz(2.3073763) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(0.63876736) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9546982) q[1];
sx q[1];
rz(-0.3317301) q[1];
sx q[1];
rz(-1.3518672) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7616368) q[3];
sx q[3];
rz(-1.8075426) q[3];
sx q[3];
rz(-0.56425795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80031359) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(-3.0467765) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.512758) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(-1.7839292) q[2];
sx q[2];
rz(-1.6518946) q[2];
sx q[2];
rz(1.2988731) q[2];
rz(2.3306866) q[3];
sx q[3];
rz(-1.9716284) q[3];
sx q[3];
rz(1.6104094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];