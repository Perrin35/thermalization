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
rz(-2.0353844) q[0];
sx q[0];
rz(2.455403) q[0];
sx q[0];
rz(10.578293) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(-2.6931813) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7767186) q[0];
sx q[0];
rz(-1.5176511) q[0];
sx q[0];
rz(2.4928635) q[0];
rz(-pi) q[1];
rz(-2.7712819) q[2];
sx q[2];
rz(-2.8680621) q[2];
sx q[2];
rz(0.47765484) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55576217) q[1];
sx q[1];
rz(-0.96802789) q[1];
sx q[1];
rz(0.65048154) q[1];
rz(-1.6187864) q[3];
sx q[3];
rz(-2.0715044) q[3];
sx q[3];
rz(-1.2402616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51556921) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(0.62421978) q[2];
rz(-0.37912399) q[3];
sx q[3];
rz(-2.1513394) q[3];
sx q[3];
rz(-0.24851255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9666331) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(2.1061184) q[0];
rz(-1.8677208) q[1];
sx q[1];
rz(-2.404411) q[1];
sx q[1];
rz(0.48286352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.14014) q[0];
sx q[0];
rz(-1.4743544) q[0];
sx q[0];
rz(1.6167377) q[0];
x q[1];
rz(0.16629433) q[2];
sx q[2];
rz(-1.6876843) q[2];
sx q[2];
rz(2.6720195) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4299996) q[1];
sx q[1];
rz(-1.5614127) q[1];
sx q[1];
rz(1.1389772) q[1];
rz(2.9008174) q[3];
sx q[3];
rz(-2.8117883) q[3];
sx q[3];
rz(-2.2036116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.024293385) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(1.8710322) q[2];
rz(-0.67000669) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(0.89699927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(0.73905149) q[0];
rz(0.96145472) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(2.3846073) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3630545) q[0];
sx q[0];
rz(-0.86150384) q[0];
sx q[0];
rz(0.34687931) q[0];
x q[1];
rz(-2.6967825) q[2];
sx q[2];
rz(-2.399657) q[2];
sx q[2];
rz(-0.43659376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0464335) q[1];
sx q[1];
rz(-0.94645247) q[1];
sx q[1];
rz(-1.4051564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7875761) q[3];
sx q[3];
rz(-2.8269973) q[3];
sx q[3];
rz(-1.5217239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6387393) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(2.4364831) q[2];
rz(2.8201568) q[3];
sx q[3];
rz(-1.0175846) q[3];
sx q[3];
rz(-2.7307935) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0099156378) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(-1.9158844) q[0];
rz(1.2340087) q[1];
sx q[1];
rz(-1.0154513) q[1];
sx q[1];
rz(3.0753678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8787251) q[0];
sx q[0];
rz(-1.6604794) q[0];
sx q[0];
rz(1.32442) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29860626) q[2];
sx q[2];
rz(-2.3176607) q[2];
sx q[2];
rz(-0.9048942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5410616) q[1];
sx q[1];
rz(-2.5503073) q[1];
sx q[1];
rz(2.1670684) q[1];
x q[2];
rz(-2.7755967) q[3];
sx q[3];
rz(-1.4533236) q[3];
sx q[3];
rz(0.42607301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0464728) q[2];
sx q[2];
rz(-2.1472223) q[2];
sx q[2];
rz(0.44060102) q[2];
rz(2.1353841) q[3];
sx q[3];
rz(-1.3738084) q[3];
sx q[3];
rz(-0.83427507) q[3];
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
rz(0.42136583) q[0];
sx q[0];
rz(-2.328673) q[0];
sx q[0];
rz(-1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.2200049) q[1];
sx q[1];
rz(1.6857326) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2746873) q[0];
sx q[0];
rz(-2.2969249) q[0];
sx q[0];
rz(-0.19498904) q[0];
rz(-1.6873932) q[2];
sx q[2];
rz(-1.3449838) q[2];
sx q[2];
rz(-1.9436398) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29042915) q[1];
sx q[1];
rz(-2.4896087) q[1];
sx q[1];
rz(3.1080721) q[1];
rz(-pi) q[2];
rz(-1.7248575) q[3];
sx q[3];
rz(-0.96821456) q[3];
sx q[3];
rz(-1.3101206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27318925) q[2];
sx q[2];
rz(-1.5634147) q[2];
sx q[2];
rz(-2.2301162) q[2];
rz(-0.8738001) q[3];
sx q[3];
rz(-2.3995212) q[3];
sx q[3];
rz(-3.1349283) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28354302) q[0];
sx q[0];
rz(-0.25074211) q[0];
sx q[0];
rz(-0.13033303) q[0];
rz(-2.6538972) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(-0.74388751) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34602964) q[0];
sx q[0];
rz(-1.6206121) q[0];
sx q[0];
rz(-1.6532142) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1962847) q[2];
sx q[2];
rz(-1.4577366) q[2];
sx q[2];
rz(0.94719749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53092693) q[1];
sx q[1];
rz(-1.0333038) q[1];
sx q[1];
rz(2.8100138) q[1];
rz(-pi) q[2];
rz(-1.7508932) q[3];
sx q[3];
rz(-1.0922179) q[3];
sx q[3];
rz(-0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0222212) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(-2.1807097) q[2];
rz(0.39572257) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(3.0254288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457086) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(1.047026) q[0];
rz(-2.537435) q[1];
sx q[1];
rz(-1.2774757) q[1];
sx q[1];
rz(-0.39271694) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0020802845) q[0];
sx q[0];
rz(-0.89186984) q[0];
sx q[0];
rz(-0.2261136) q[0];
rz(0.95177841) q[2];
sx q[2];
rz(-1.5238395) q[2];
sx q[2];
rz(-2.6359735) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8046569) q[1];
sx q[1];
rz(-1.5976904) q[1];
sx q[1];
rz(-2.8662445) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32870558) q[3];
sx q[3];
rz(-1.3188601) q[3];
sx q[3];
rz(-3.1303034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.087223209) q[2];
sx q[2];
rz(-1.1944218) q[2];
sx q[2];
rz(1.3090022) q[2];
rz(1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(-0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.29276174) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(0.30717474) q[0];
rz(1.81709) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(-0.90075341) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15316111) q[0];
sx q[0];
rz(-1.5539196) q[0];
sx q[0];
rz(2.8794226) q[0];
x q[1];
rz(-1.0020761) q[2];
sx q[2];
rz(-0.7373215) q[2];
sx q[2];
rz(-1.3739283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0010502) q[1];
sx q[1];
rz(-1.2549572) q[1];
sx q[1];
rz(-1.9204101) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7124655) q[3];
sx q[3];
rz(-0.9441388) q[3];
sx q[3];
rz(0.10654813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8073392) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(-0.61484289) q[2];
rz(-2.0452512) q[3];
sx q[3];
rz(-2.5494826) q[3];
sx q[3];
rz(2.443327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7480302) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(-0.50931859) q[0];
rz(-0.18109426) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(0.96955713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80705331) q[0];
sx q[0];
rz(-1.7248123) q[0];
sx q[0];
rz(-0.75169433) q[0];
rz(3.0907057) q[2];
sx q[2];
rz(-2.0774088) q[2];
sx q[2];
rz(-1.2801054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4050715) q[1];
sx q[1];
rz(-2.7847291) q[1];
sx q[1];
rz(2.8282911) q[1];
rz(-1.182231) q[3];
sx q[3];
rz(-2.3747184) q[3];
sx q[3];
rz(-2.9076613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2481044) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-2.060037) q[2];
rz(-0.23165101) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(-2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3807826) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(2.494452) q[0];
rz(-2.6864247) q[1];
sx q[1];
rz(-1.7166694) q[1];
sx q[1];
rz(0.761935) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1200346) q[0];
sx q[0];
rz(-1.2823449) q[0];
sx q[0];
rz(-0.7201654) q[0];
rz(2.5013431) q[2];
sx q[2];
rz(-0.92776042) q[2];
sx q[2];
rz(-1.6563479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.476452) q[1];
sx q[1];
rz(-2.2287205) q[1];
sx q[1];
rz(1.5275147) q[1];
x q[2];
rz(1.1033789) q[3];
sx q[3];
rz(-1.5639135) q[3];
sx q[3];
rz(2.9228743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.053293856) q[2];
sx q[2];
rz(-0.21685728) q[2];
sx q[2];
rz(-2.2376412) q[2];
rz(1.5229185) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91697964) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(-3.0147973) q[1];
sx q[1];
rz(-2.3190111) q[1];
sx q[1];
rz(0.93228985) q[1];
rz(-2.8895072) q[2];
sx q[2];
rz(-2.0093976) q[2];
sx q[2];
rz(-1.2653399) q[2];
rz(2.7115962) q[3];
sx q[3];
rz(-1.5734133) q[3];
sx q[3];
rz(-1.6071241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
