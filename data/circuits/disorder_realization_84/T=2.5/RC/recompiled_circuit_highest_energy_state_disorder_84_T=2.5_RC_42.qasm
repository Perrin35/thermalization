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
rz(2.5166729) q[0];
sx q[0];
rz(-1.807037) q[0];
sx q[0];
rz(-3.0883489) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(6.1558131) q[1];
sx q[1];
rz(11.131412) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2086528) q[0];
sx q[0];
rz(-2.8253252) q[0];
sx q[0];
rz(0.44751944) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26726079) q[2];
sx q[2];
rz(-1.6858941) q[2];
sx q[2];
rz(-0.54265825) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8813144) q[1];
sx q[1];
rz(-1.5173755) q[1];
sx q[1];
rz(0.39398663) q[1];
rz(0.93115689) q[3];
sx q[3];
rz(-1.1619029) q[3];
sx q[3];
rz(-2.3774022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3367553) q[2];
sx q[2];
rz(-1.7982322) q[2];
sx q[2];
rz(-2.8678144) q[2];
rz(-1.2182419) q[3];
sx q[3];
rz(-0.67353407) q[3];
sx q[3];
rz(-2.8555253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.022920595) q[0];
sx q[0];
rz(-0.76347041) q[0];
sx q[0];
rz(0.77962312) q[0];
rz(-2.4769705) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(1.2453311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0459291) q[0];
sx q[0];
rz(-1.6324658) q[0];
sx q[0];
rz(-1.685919) q[0];
rz(1.677956) q[2];
sx q[2];
rz(-0.71427155) q[2];
sx q[2];
rz(2.799054) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7088931) q[1];
sx q[1];
rz(-1.3586167) q[1];
sx q[1];
rz(0.062122155) q[1];
rz(-pi) q[2];
rz(2.1297203) q[3];
sx q[3];
rz(-1.3959612) q[3];
sx q[3];
rz(1.2993985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43210426) q[2];
sx q[2];
rz(-0.69807845) q[2];
sx q[2];
rz(-3.091605) q[2];
rz(-1.4738119) q[3];
sx q[3];
rz(-1.7260909) q[3];
sx q[3];
rz(0.93562359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1525928) q[0];
sx q[0];
rz(-0.86243668) q[0];
sx q[0];
rz(-1.1784026) q[0];
rz(-1.1475457) q[1];
sx q[1];
rz(-2.9541364) q[1];
sx q[1];
rz(0.60290927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99166223) q[0];
sx q[0];
rz(-1.7310744) q[0];
sx q[0];
rz(0.10326881) q[0];
x q[1];
rz(-0.96430594) q[2];
sx q[2];
rz(-0.78769257) q[2];
sx q[2];
rz(0.76157969) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29290043) q[1];
sx q[1];
rz(-2.22977) q[1];
sx q[1];
rz(2.1007358) q[1];
rz(-pi) q[2];
rz(-2.1985377) q[3];
sx q[3];
rz(-2.3585417) q[3];
sx q[3];
rz(2.8179226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9511562) q[2];
sx q[2];
rz(-0.78593212) q[2];
sx q[2];
rz(2.7583165) q[2];
rz(-1.8303309) q[3];
sx q[3];
rz(-2.0537328) q[3];
sx q[3];
rz(0.12152984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.424778) q[0];
sx q[0];
rz(-2.1489693) q[0];
sx q[0];
rz(0.92856652) q[0];
rz(0.83944744) q[1];
sx q[1];
rz(-1.8239559) q[1];
sx q[1];
rz(-2.3323257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033484785) q[0];
sx q[0];
rz(-1.9376457) q[0];
sx q[0];
rz(-0.97292329) q[0];
rz(-pi) q[1];
rz(-2.7847084) q[2];
sx q[2];
rz(-0.45512558) q[2];
sx q[2];
rz(-2.9659716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2553808) q[1];
sx q[1];
rz(-1.501923) q[1];
sx q[1];
rz(-1.3248454) q[1];
rz(-0.92692848) q[3];
sx q[3];
rz(-1.2371105) q[3];
sx q[3];
rz(2.9179171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75260085) q[2];
sx q[2];
rz(-1.6099124) q[2];
sx q[2];
rz(1.4468225) q[2];
rz(-1.1270479) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(2.5126422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.9160354) q[0];
sx q[0];
rz(-1.9434403) q[0];
sx q[0];
rz(1.7115364) q[0];
rz(-2.9043708) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(-0.29475862) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4859516) q[0];
sx q[0];
rz(-1.7136586) q[0];
sx q[0];
rz(2.7014707) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7698603) q[2];
sx q[2];
rz(-0.9424389) q[2];
sx q[2];
rz(-0.5231176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.393906) q[1];
sx q[1];
rz(-1.1562041) q[1];
sx q[1];
rz(1.0769597) q[1];
x q[2];
rz(1.7196017) q[3];
sx q[3];
rz(-2.1396881) q[3];
sx q[3];
rz(2.1619145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1416867) q[2];
sx q[2];
rz(-0.20076951) q[2];
sx q[2];
rz(-3.0086009) q[2];
rz(-0.76987949) q[3];
sx q[3];
rz(-1.2917638) q[3];
sx q[3];
rz(-2.8225115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.28984508) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(-1.7293365) q[0];
rz(-0.051636592) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(-2.9488865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7190921) q[0];
sx q[0];
rz(-1.0570881) q[0];
sx q[0];
rz(-1.5925608) q[0];
rz(0.66484837) q[2];
sx q[2];
rz(-2.2908205) q[2];
sx q[2];
rz(0.31535044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51260692) q[1];
sx q[1];
rz(-0.64388212) q[1];
sx q[1];
rz(-0.75276796) q[1];
x q[2];
rz(-0.64151836) q[3];
sx q[3];
rz(-1.8810442) q[3];
sx q[3];
rz(1.9812685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.073033832) q[2];
sx q[2];
rz(-1.8517588) q[2];
sx q[2];
rz(-1.3812836) q[2];
rz(-0.51501385) q[3];
sx q[3];
rz(-2.2541855) q[3];
sx q[3];
rz(1.7611586) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16159049) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(2.7951151) q[0];
rz(2.1222291) q[1];
sx q[1];
rz(-2.4788224) q[1];
sx q[1];
rz(1.4541218) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0270868) q[0];
sx q[0];
rz(-2.1505304) q[0];
sx q[0];
rz(1.0282577) q[0];
rz(2.2367661) q[2];
sx q[2];
rz(-0.77205333) q[2];
sx q[2];
rz(-1.1386629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4291961) q[1];
sx q[1];
rz(-2.0919261) q[1];
sx q[1];
rz(-2.8217535) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4164657) q[3];
sx q[3];
rz(-1.6701506) q[3];
sx q[3];
rz(-1.1442483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5399897) q[2];
sx q[2];
rz(-1.8113965) q[2];
sx q[2];
rz(2.9023564) q[2];
rz(2.7573977) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(2.189883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20188986) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(1.3772759) q[0];
rz(-1.2205623) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(-2.9753704) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3275422) q[0];
sx q[0];
rz(-0.53314994) q[0];
sx q[0];
rz(2.0943805) q[0];
rz(1.8883838) q[2];
sx q[2];
rz(-0.86195213) q[2];
sx q[2];
rz(2.6203757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5095013) q[1];
sx q[1];
rz(-2.9281514) q[1];
sx q[1];
rz(-0.14519338) q[1];
x q[2];
rz(1.1508302) q[3];
sx q[3];
rz(-1.3356497) q[3];
sx q[3];
rz(1.8311892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.347747) q[2];
sx q[2];
rz(-2.2913427) q[2];
sx q[2];
rz(1.0469077) q[2];
rz(1.8868123) q[3];
sx q[3];
rz(-1.0898277) q[3];
sx q[3];
rz(1.8449529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60085249) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(-1.1131713) q[0];
rz(-0.81740776) q[1];
sx q[1];
rz(-1.6956885) q[1];
sx q[1];
rz(2.7730952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2558852) q[0];
sx q[0];
rz(-1.8428546) q[0];
sx q[0];
rz(1.9291124) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13780528) q[2];
sx q[2];
rz(-0.74713444) q[2];
sx q[2];
rz(-1.8835619) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.52942136) q[1];
sx q[1];
rz(-2.2670806) q[1];
sx q[1];
rz(1.2207635) q[1];
x q[2];
rz(2.1300674) q[3];
sx q[3];
rz(-0.91010053) q[3];
sx q[3];
rz(-1.6370156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9937399) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.7842133) q[2];
rz(-2.3944858) q[3];
sx q[3];
rz(-2.1785469) q[3];
sx q[3];
rz(-0.67556206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6152182) q[0];
sx q[0];
rz(-2.6850061) q[0];
sx q[0];
rz(-1.0427465) q[0];
rz(-2.3710947) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(-0.76046336) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2337991) q[0];
sx q[0];
rz(-1.9007508) q[0];
sx q[0];
rz(2.8889662) q[0];
rz(-1.3907688) q[2];
sx q[2];
rz(-1.5329201) q[2];
sx q[2];
rz(-1.2779209) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9925966) q[1];
sx q[1];
rz(-1.9933874) q[1];
sx q[1];
rz(1.6572957) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46871878) q[3];
sx q[3];
rz(-0.7112452) q[3];
sx q[3];
rz(-0.48392344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.20327917) q[2];
sx q[2];
rz(-0.61050296) q[2];
sx q[2];
rz(0.40943134) q[2];
rz(1.5583386) q[3];
sx q[3];
rz(-0.46054545) q[3];
sx q[3];
rz(2.515365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4938477) q[0];
sx q[0];
rz(-1.0418325) q[0];
sx q[0];
rz(-2.7000725) q[0];
rz(1.582675) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(-2.6263993) q[2];
sx q[2];
rz(-0.11820785) q[2];
sx q[2];
rz(-1.6109656) q[2];
rz(1.7033475) q[3];
sx q[3];
rz(-1.3302531) q[3];
sx q[3];
rz(2.1147685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
