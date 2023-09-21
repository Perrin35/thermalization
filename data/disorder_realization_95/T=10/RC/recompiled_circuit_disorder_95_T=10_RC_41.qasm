OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(-2.2440417) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25047725) q[0];
sx q[0];
rz(-2.5076137) q[0];
sx q[0];
rz(-1.4527713) q[0];
rz(1.1041553) q[2];
sx q[2];
rz(-0.68744126) q[2];
sx q[2];
rz(-1.825037) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0734288) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(2.8627002) q[1];
rz(-pi) q[2];
rz(2.3834121) q[3];
sx q[3];
rz(-1.7054134) q[3];
sx q[3];
rz(-1.7077703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.477318) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(-2.409639) q[2];
rz(-2.1814363) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(2.2629471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4449578) q[0];
sx q[0];
rz(-0.96224552) q[0];
sx q[0];
rz(-0.47136013) q[0];
x q[1];
rz(-0.61848817) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(-2.6478812) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6138184) q[1];
sx q[1];
rz(-1.4973745) q[1];
sx q[1];
rz(-0.41927494) q[1];
x q[2];
rz(-3.0733172) q[3];
sx q[3];
rz(-1.2618511) q[3];
sx q[3];
rz(0.49648778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(-2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-1.1211959) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58223984) q[0];
sx q[0];
rz(-1.9358578) q[0];
sx q[0];
rz(2.3854726) q[0];
rz(-pi) q[1];
rz(-1.7240702) q[2];
sx q[2];
rz(-2.7445265) q[2];
sx q[2];
rz(3.0561471) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0507625) q[1];
sx q[1];
rz(-1.275626) q[1];
sx q[1];
rz(-0.34543085) q[1];
rz(-0.49384533) q[3];
sx q[3];
rz(-1.4673127) q[3];
sx q[3];
rz(0.34085694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(-0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-2.3535368) q[0];
rz(-2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(3.025211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6560293) q[0];
sx q[0];
rz(-2.2040327) q[0];
sx q[0];
rz(-1.602519) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38168455) q[2];
sx q[2];
rz(-1.0978062) q[2];
sx q[2];
rz(0.7009398) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9434005) q[1];
sx q[1];
rz(-1.0775837) q[1];
sx q[1];
rz(-2.0342159) q[1];
rz(3.0451123) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(-1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(0.3616412) q[3];
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
rz(pi/2) q[3];
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
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.9929569) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900433) q[0];
sx q[0];
rz(-0.94928375) q[0];
sx q[0];
rz(-3.050699) q[0];
x q[1];
rz(-2.6905641) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(0.23652467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3061805) q[1];
sx q[1];
rz(-0.40459834) q[1];
sx q[1];
rz(-0.49680357) q[1];
rz(-pi) q[2];
rz(-0.72539056) q[3];
sx q[3];
rz(-1.8728421) q[3];
sx q[3];
rz(0.52291742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1308412) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(-2.1095236) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-2.9352303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9128742) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(2.3939783) q[0];
rz(2.3328853) q[2];
sx q[2];
rz(-0.84918298) q[2];
sx q[2];
rz(-1.2869814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3316162) q[1];
sx q[1];
rz(-2.9990701) q[1];
sx q[1];
rz(-2.317821) q[1];
rz(1.0422802) q[3];
sx q[3];
rz(-0.37502608) q[3];
sx q[3];
rz(-3.022775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(-2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966184) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.2111838) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1887218) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(-2.1923724) q[0];
rz(2.6518875) q[2];
sx q[2];
rz(-1.6929132) q[2];
sx q[2];
rz(0.44039886) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3019575) q[1];
sx q[1];
rz(-0.17050276) q[1];
sx q[1];
rz(-1.9819928) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1189506) q[3];
sx q[3];
rz(-0.65330905) q[3];
sx q[3];
rz(2.2484231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3367735) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(-3.1170735) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0054758469) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(2.1210282) q[0];
rz(2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.221009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2894665) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(0.37129398) q[0];
rz(2.3983725) q[2];
sx q[2];
rz(-0.23822242) q[2];
sx q[2];
rz(1.5526349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.358566) q[1];
sx q[1];
rz(-1.6811803) q[1];
sx q[1];
rz(-0.13659887) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2113167) q[3];
sx q[3];
rz(-0.7822789) q[3];
sx q[3];
rz(0.6560916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(-2.4364046) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(1.5195297) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(0.12577122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3868956) q[0];
sx q[0];
rz(-1.3820952) q[0];
sx q[0];
rz(2.9625091) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52351064) q[2];
sx q[2];
rz(-2.7611809) q[2];
sx q[2];
rz(-1.3695804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51024918) q[1];
sx q[1];
rz(-1.2878839) q[1];
sx q[1];
rz(0.81088539) q[1];
x q[2];
rz(-2.845876) q[3];
sx q[3];
rz(-1.3774646) q[3];
sx q[3];
rz(2.8545692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(-2.4323145) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-0.15670776) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9897292) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(-1.1402003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4338019) q[0];
sx q[0];
rz(-1.3085438) q[0];
sx q[0];
rz(0.10284013) q[0];
rz(-0.26009772) q[2];
sx q[2];
rz(-1.1220699) q[2];
sx q[2];
rz(-1.5990433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4789341) q[1];
sx q[1];
rz(-1.7591898) q[1];
sx q[1];
rz(-1.6608095) q[1];
rz(1.8885918) q[3];
sx q[3];
rz(-2.3386049) q[3];
sx q[3];
rz(-3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(0.89861384) q[2];
rz(3.1344154) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(2.6208411) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(-2.043914) q[2];
sx q[2];
rz(-2.0156779) q[2];
sx q[2];
rz(1.426522) q[2];
rz(0.58581523) q[3];
sx q[3];
rz(-1.2231493) q[3];
sx q[3];
rz(2.3412658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];