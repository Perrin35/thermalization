OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(0.97283483) q[1];
sx q[1];
rz(-1.4714779) q[1];
sx q[1];
rz(0.29247984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7115294) q[0];
sx q[0];
rz(-1.1769413) q[0];
sx q[0];
rz(0.060056134) q[0];
rz(-pi) q[1];
rz(1.6220665) q[2];
sx q[2];
rz(-0.79203696) q[2];
sx q[2];
rz(-2.771332) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3136183) q[1];
sx q[1];
rz(-2.5458114) q[1];
sx q[1];
rz(1.4890563) q[1];
rz(-pi) q[2];
rz(1.5745893) q[3];
sx q[3];
rz(-2.6677701) q[3];
sx q[3];
rz(-0.68236085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4108489) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(2.6921819) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(-1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(0.74203062) q[0];
rz(-1.6702601) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.7785633) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3516386) q[0];
sx q[0];
rz(-2.3839256) q[0];
sx q[0];
rz(0.58654465) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6724733) q[2];
sx q[2];
rz(-1.763478) q[2];
sx q[2];
rz(-2.0872781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7098397) q[1];
sx q[1];
rz(-1.5515986) q[1];
sx q[1];
rz(-0.33961105) q[1];
rz(-pi) q[2];
rz(-0.95008534) q[3];
sx q[3];
rz(-2.8182497) q[3];
sx q[3];
rz(-2.8107373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8403975) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.8691241) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0740046) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(-2.8787676) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6806157) q[0];
sx q[0];
rz(-1.5605698) q[0];
sx q[0];
rz(0.015120487) q[0];
rz(-pi) q[1];
rz(-1.8112196) q[2];
sx q[2];
rz(-1.5449107) q[2];
sx q[2];
rz(2.7607069) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0318109) q[1];
sx q[1];
rz(-0.33966741) q[1];
sx q[1];
rz(-0.32082816) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34706195) q[3];
sx q[3];
rz(-2.2357781) q[3];
sx q[3];
rz(-2.5393328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1245023) q[2];
sx q[2];
rz(-1.7586781) q[2];
sx q[2];
rz(-0.18033218) q[2];
rz(-2.5056433) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(-0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.3746049) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(2.1069353) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67544115) q[0];
sx q[0];
rz(-3.0794641) q[0];
sx q[0];
rz(0.32390578) q[0];
rz(-pi) q[1];
rz(1.9070508) q[2];
sx q[2];
rz(-0.14306919) q[2];
sx q[2];
rz(-2.2640995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82995854) q[1];
sx q[1];
rz(-1.2882075) q[1];
sx q[1];
rz(2.4294873) q[1];
x q[2];
rz(-0.33200522) q[3];
sx q[3];
rz(-2.6453291) q[3];
sx q[3];
rz(0.23976025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.93306142) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(2.9053524) q[2];
rz(1.9832206) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(-0.90233666) q[0];
rz(2.2205655) q[1];
sx q[1];
rz(-0.61704707) q[1];
sx q[1];
rz(-2.5193118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392004) q[0];
sx q[0];
rz(-1.3883739) q[0];
sx q[0];
rz(1.8586041) q[0];
x q[1];
rz(1.4786517) q[2];
sx q[2];
rz(-0.92458692) q[2];
sx q[2];
rz(0.19942936) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29187782) q[1];
sx q[1];
rz(-1.7397225) q[1];
sx q[1];
rz(3.0506163) q[1];
rz(-1.9641818) q[3];
sx q[3];
rz(-1.1227566) q[3];
sx q[3];
rz(-0.084669948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(1.4422669) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(3.0969627) q[0];
rz(1.0021098) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(-0.034428509) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32193004) q[0];
sx q[0];
rz(-2.9557883) q[0];
sx q[0];
rz(-1.6447322) q[0];
rz(-3.0197633) q[2];
sx q[2];
rz(-2.0707154) q[2];
sx q[2];
rz(-0.027539754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5163102) q[1];
sx q[1];
rz(-2.1045661) q[1];
sx q[1];
rz(2.7620402) q[1];
rz(-pi) q[2];
rz(0.70528443) q[3];
sx q[3];
rz(-1.4449638) q[3];
sx q[3];
rz(-0.82334405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3520711) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(-2.1748523) q[2];
rz(-3.1294075) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(-0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352585) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(-0.99408856) q[0];
rz(3.0864691) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(-2.4023043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2544884) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(-1.7083005) q[0];
rz(-pi) q[1];
rz(-2.0544858) q[2];
sx q[2];
rz(-2.0927883) q[2];
sx q[2];
rz(2.1453478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48078295) q[1];
sx q[1];
rz(-1.0479095) q[1];
sx q[1];
rz(-1.3516264) q[1];
x q[2];
rz(2.5218511) q[3];
sx q[3];
rz(-1.8078139) q[3];
sx q[3];
rz(1.2902416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5560975) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(2.901315) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69650841) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(0.58404303) q[0];
rz(0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(-2.8631794) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42785545) q[0];
sx q[0];
rz(-0.58361485) q[0];
sx q[0];
rz(1.3419271) q[0];
x q[1];
rz(1.3827219) q[2];
sx q[2];
rz(-1.1838059) q[2];
sx q[2];
rz(2.219971) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1380996) q[1];
sx q[1];
rz(-1.3521191) q[1];
sx q[1];
rz(-1.9298645) q[1];
rz(-pi) q[2];
rz(-1.9592459) q[3];
sx q[3];
rz(-2.5591345) q[3];
sx q[3];
rz(-2.9782481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(-0.38796866) q[2];
rz(1.3502454) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(-0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(0.073721185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6442935) q[0];
sx q[0];
rz(-1.6719581) q[0];
sx q[0];
rz(-1.1343207) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6905626) q[2];
sx q[2];
rz(-1.6511859) q[2];
sx q[2];
rz(0.62270852) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24450609) q[1];
sx q[1];
rz(-0.58468854) q[1];
sx q[1];
rz(-0.067824407) q[1];
rz(-0.1707465) q[3];
sx q[3];
rz(-2.5617544) q[3];
sx q[3];
rz(0.29614007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5471197) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(2.7129042) q[2];
rz(-1.3171014) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(-1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4093032) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(1.0428585) q[0];
rz(-1.6304784) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(-1.7932549) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6207328) q[0];
sx q[0];
rz(-0.11611406) q[0];
sx q[0];
rz(1.3470879) q[0];
rz(-pi) q[1];
rz(1.6595608) q[2];
sx q[2];
rz(-2.4030493) q[2];
sx q[2];
rz(0.86631394) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18689449) q[1];
sx q[1];
rz(-0.3317301) q[1];
sx q[1];
rz(1.7897254) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24095778) q[3];
sx q[3];
rz(-1.3853419) q[3];
sx q[3];
rz(1.0518187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(3.0467765) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(-0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.6288347) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(-1.6013153) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(-1.3576635) q[2];
sx q[2];
rz(-1.4896981) q[2];
sx q[2];
rz(-1.8427195) q[2];
rz(0.52900984) q[3];
sx q[3];
rz(-0.88376868) q[3];
sx q[3];
rz(0.39467011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
