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
rz(1.6871356) q[0];
sx q[0];
rz(-1.1633101) q[0];
sx q[0];
rz(1.8736725) q[0];
rz(1.6549702) q[1];
sx q[1];
rz(4.4284664) q[1];
sx q[1];
rz(9.5944302) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96972966) q[0];
sx q[0];
rz(-0.61530441) q[0];
sx q[0];
rz(0.91795532) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1450139) q[2];
sx q[2];
rz(-0.55934956) q[2];
sx q[2];
rz(2.1970791) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7387661) q[1];
sx q[1];
rz(-1.9662274) q[1];
sx q[1];
rz(-0.02301245) q[1];
x q[2];
rz(0.98244169) q[3];
sx q[3];
rz(-0.37068493) q[3];
sx q[3];
rz(-1.4511758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22592813) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(-0.020326745) q[2];
rz(-2.9018371) q[3];
sx q[3];
rz(-0.25529796) q[3];
sx q[3];
rz(0.83893004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13424419) q[0];
sx q[0];
rz(-0.61779314) q[0];
sx q[0];
rz(-1.0963305) q[0];
rz(-1.4758551) q[1];
sx q[1];
rz(-1.9225537) q[1];
sx q[1];
rz(2.347167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5223094) q[0];
sx q[0];
rz(-1.8849775) q[0];
sx q[0];
rz(0.35616014) q[0];
x q[1];
rz(-2.3874823) q[2];
sx q[2];
rz(-1.1271141) q[2];
sx q[2];
rz(0.59690969) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5880044) q[1];
sx q[1];
rz(-0.4376227) q[1];
sx q[1];
rz(-0.2511843) q[1];
rz(-pi) q[2];
rz(-0.81476029) q[3];
sx q[3];
rz(-1.4346806) q[3];
sx q[3];
rz(2.5825115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0558001) q[2];
sx q[2];
rz(-1.5726568) q[2];
sx q[2];
rz(-2.9158578) q[2];
rz(0.24383946) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(-0.17914151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8826411) q[0];
sx q[0];
rz(-1.8283586) q[0];
sx q[0];
rz(0.42695811) q[0];
rz(0.34307617) q[1];
sx q[1];
rz(-1.1767358) q[1];
sx q[1];
rz(-2.8899736) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31684012) q[0];
sx q[0];
rz(-1.7946294) q[0];
sx q[0];
rz(0.59455183) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.919056) q[2];
sx q[2];
rz(-2.3250569) q[2];
sx q[2];
rz(2.726462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15089825) q[1];
sx q[1];
rz(-0.37726918) q[1];
sx q[1];
rz(2.6713085) q[1];
x q[2];
rz(0.86029025) q[3];
sx q[3];
rz(-1.6087469) q[3];
sx q[3];
rz(0.52805985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.130927) q[2];
sx q[2];
rz(-1.6341354) q[2];
sx q[2];
rz(-0.44535401) q[2];
rz(1.2238067) q[3];
sx q[3];
rz(-1.8755951) q[3];
sx q[3];
rz(0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39684108) q[0];
sx q[0];
rz(-0.80533177) q[0];
sx q[0];
rz(1.1591563) q[0];
rz(-1.9388439) q[1];
sx q[1];
rz(-1.9614599) q[1];
sx q[1];
rz(0.73701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4241959) q[0];
sx q[0];
rz(-1.5079466) q[0];
sx q[0];
rz(2.8934188) q[0];
rz(-pi) q[1];
rz(-1.1221755) q[2];
sx q[2];
rz(-2.4146842) q[2];
sx q[2];
rz(2.2688933) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9056978) q[1];
sx q[1];
rz(-1.7445843) q[1];
sx q[1];
rz(2.5008965) q[1];
x q[2];
rz(2.9748671) q[3];
sx q[3];
rz(-1.5136216) q[3];
sx q[3];
rz(-1.0852607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1272588) q[2];
sx q[2];
rz(-2.2517683) q[2];
sx q[2];
rz(-0.8117525) q[2];
rz(0.74770606) q[3];
sx q[3];
rz(-0.85993189) q[3];
sx q[3];
rz(0.060613304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6478445) q[0];
sx q[0];
rz(-1.0890549) q[0];
sx q[0];
rz(-1.7895948) q[0];
rz(-0.79212517) q[1];
sx q[1];
rz(-0.89901662) q[1];
sx q[1];
rz(-2.1314714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5631752) q[0];
sx q[0];
rz(-1.0437328) q[0];
sx q[0];
rz(0.38946797) q[0];
rz(-pi) q[1];
rz(0.13212684) q[2];
sx q[2];
rz(-0.39574049) q[2];
sx q[2];
rz(-2.6504315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5252851) q[1];
sx q[1];
rz(-1.0152363) q[1];
sx q[1];
rz(1.4772593) q[1];
rz(-pi) q[2];
rz(-1.7351446) q[3];
sx q[3];
rz(-1.4557299) q[3];
sx q[3];
rz(0.0054520741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1048364) q[2];
sx q[2];
rz(-0.75883055) q[2];
sx q[2];
rz(1.7581615) q[2];
rz(-3.1116327) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(1.8647319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.5112011) q[0];
sx q[0];
rz(-0.41655219) q[0];
sx q[0];
rz(-3.0101486) q[0];
rz(0.99880544) q[1];
sx q[1];
rz(-1.1591594) q[1];
sx q[1];
rz(1.3712032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72238647) q[0];
sx q[0];
rz(-2.6372069) q[0];
sx q[0];
rz(0.76865159) q[0];
rz(-0.24662416) q[2];
sx q[2];
rz(-0.70136753) q[2];
sx q[2];
rz(-1.3878617) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35624337) q[1];
sx q[1];
rz(-1.8438854) q[1];
sx q[1];
rz(2.9778984) q[1];
rz(-pi) q[2];
rz(-0.71089069) q[3];
sx q[3];
rz(-1.5734473) q[3];
sx q[3];
rz(0.34429911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3580631) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(2.0520468) q[2];
rz(0.36695668) q[3];
sx q[3];
rz(-0.4042545) q[3];
sx q[3];
rz(0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628196) q[0];
sx q[0];
rz(-2.3371526) q[0];
sx q[0];
rz(0.87271571) q[0];
rz(-1.577042) q[1];
sx q[1];
rz(-1.6460452) q[1];
sx q[1];
rz(2.4094792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4622065) q[0];
sx q[0];
rz(-2.4024978) q[0];
sx q[0];
rz(-2.8406891) q[0];
rz(-2.9896834) q[2];
sx q[2];
rz(-0.6912125) q[2];
sx q[2];
rz(0.37746261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9804025) q[1];
sx q[1];
rz(-1.727961) q[1];
sx q[1];
rz(2.0445301) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0454265) q[3];
sx q[3];
rz(-1.3447666) q[3];
sx q[3];
rz(-1.9440252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6096036) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(-0.58094376) q[2];
rz(-2.0161435) q[3];
sx q[3];
rz(-1.9476798) q[3];
sx q[3];
rz(1.1552936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67227143) q[0];
sx q[0];
rz(-0.10675616) q[0];
sx q[0];
rz(-2.971055) q[0];
rz(1.8960309) q[1];
sx q[1];
rz(-1.6163369) q[1];
sx q[1];
rz(2.4571498) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67407284) q[0];
sx q[0];
rz(-1.0985803) q[0];
sx q[0];
rz(-1.1961106) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6587517) q[2];
sx q[2];
rz(-1.4582658) q[2];
sx q[2];
rz(-1.4275328) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0674921) q[1];
sx q[1];
rz(-1.3252186) q[1];
sx q[1];
rz(-0.26439338) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.079250431) q[3];
sx q[3];
rz(-0.65046117) q[3];
sx q[3];
rz(-1.6975171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30274063) q[2];
sx q[2];
rz(-1.799823) q[2];
sx q[2];
rz(2.0212685) q[2];
rz(0.66425792) q[3];
sx q[3];
rz(-2.7393326) q[3];
sx q[3];
rz(2.6334488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6054194) q[0];
sx q[0];
rz(-0.37261951) q[0];
sx q[0];
rz(2.5033409) q[0];
rz(-0.0087180184) q[1];
sx q[1];
rz(-2.0254841) q[1];
sx q[1];
rz(-2.7353824) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531949) q[0];
sx q[0];
rz(-2.8450845) q[0];
sx q[0];
rz(-1.1738846) q[0];
rz(-pi) q[1];
rz(-0.052283124) q[2];
sx q[2];
rz(-2.0039275) q[2];
sx q[2];
rz(-1.0811102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9487755) q[1];
sx q[1];
rz(-1.4566167) q[1];
sx q[1];
rz(2.3523056) q[1];
x q[2];
rz(3.0681455) q[3];
sx q[3];
rz(-1.5961507) q[3];
sx q[3];
rz(-0.0097833477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0476039) q[2];
sx q[2];
rz(-2.0145907) q[2];
sx q[2];
rz(0.36063933) q[2];
rz(0.7274729) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(2.8235161) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3280846) q[0];
sx q[0];
rz(-0.94539517) q[0];
sx q[0];
rz(0.16950053) q[0];
rz(2.4318579) q[1];
sx q[1];
rz(-1.7252445) q[1];
sx q[1];
rz(1.08606) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5026131) q[0];
sx q[0];
rz(-1.3492246) q[0];
sx q[0];
rz(-0.53934877) q[0];
rz(-0.96699826) q[2];
sx q[2];
rz(-0.51568401) q[2];
sx q[2];
rz(-2.5623851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5063613) q[1];
sx q[1];
rz(-1.4725793) q[1];
sx q[1];
rz(1.0088831) q[1];
x q[2];
rz(0.38137718) q[3];
sx q[3];
rz(-1.3998195) q[3];
sx q[3];
rz(1.387763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78709948) q[2];
sx q[2];
rz(-3.0566065) q[2];
sx q[2];
rz(0.3698012) q[2];
rz(-1.1981111) q[3];
sx q[3];
rz(-1.5503649) q[3];
sx q[3];
rz(-2.2280391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030466) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(1.0176324) q[1];
sx q[1];
rz(-1.3980649) q[1];
sx q[1];
rz(-1.3921888) q[1];
rz(0.77059435) q[2];
sx q[2];
rz(-0.10100766) q[2];
sx q[2];
rz(-2.934692) q[2];
rz(-1.4403338) q[3];
sx q[3];
rz(-1.7335907) q[3];
sx q[3];
rz(2.8681267) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
