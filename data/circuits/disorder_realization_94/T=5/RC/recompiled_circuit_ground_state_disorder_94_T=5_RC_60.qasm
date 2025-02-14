OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2151467) q[0];
sx q[0];
rz(-2.4617221) q[0];
sx q[0];
rz(0.1006861) q[0];
rz(-2.6868532) q[1];
sx q[1];
rz(-2.4603031) q[1];
sx q[1];
rz(-2.1159621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5012623) q[0];
sx q[0];
rz(-2.9399031) q[0];
sx q[0];
rz(1.8224688) q[0];
rz(-2.8335613) q[2];
sx q[2];
rz(-1.5218166) q[2];
sx q[2];
rz(-0.9020976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1069965) q[1];
sx q[1];
rz(-1.6507848) q[1];
sx q[1];
rz(-0.75350113) q[1];
x q[2];
rz(2.493606) q[3];
sx q[3];
rz(-0.77666908) q[3];
sx q[3];
rz(0.69029337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2166298) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(1.487757) q[2];
rz(-1.9803068) q[3];
sx q[3];
rz(-2.1372644) q[3];
sx q[3];
rz(1.0950834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0733114) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(-0.43346369) q[0];
rz(-2.0794226) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(2.4210222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90209157) q[0];
sx q[0];
rz(-1.6974309) q[0];
sx q[0];
rz(0.037221639) q[0];
rz(-pi) q[1];
rz(-1.0373239) q[2];
sx q[2];
rz(-1.7852655) q[2];
sx q[2];
rz(2.7191424) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1432508) q[1];
sx q[1];
rz(-1.9228808) q[1];
sx q[1];
rz(1.4125173) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0508762) q[3];
sx q[3];
rz(-0.54819104) q[3];
sx q[3];
rz(-2.2238126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34708193) q[2];
sx q[2];
rz(-1.2257267) q[2];
sx q[2];
rz(-2.1573055) q[2];
rz(-2.2847564) q[3];
sx q[3];
rz(-2.555116) q[3];
sx q[3];
rz(-1.2543359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1089351) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(2.3987067) q[0];
rz(3.1318829) q[1];
sx q[1];
rz(-1.566889) q[1];
sx q[1];
rz(2.8579393) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8946202) q[0];
sx q[0];
rz(-0.63480154) q[0];
sx q[0];
rz(0.60946861) q[0];
x q[1];
rz(-2.9417324) q[2];
sx q[2];
rz(-1.4132573) q[2];
sx q[2];
rz(3.1197366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.328542) q[1];
sx q[1];
rz(-0.35085163) q[1];
sx q[1];
rz(1.0628797) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4693442) q[3];
sx q[3];
rz(-1.5939004) q[3];
sx q[3];
rz(0.35233179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0297086) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(0.65579826) q[2];
rz(2.9295975) q[3];
sx q[3];
rz(-2.5097804) q[3];
sx q[3];
rz(-1.9688212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028932171) q[0];
sx q[0];
rz(-2.6472968) q[0];
sx q[0];
rz(-1.5546881) q[0];
rz(2.9117865) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(1.302964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4939319) q[0];
sx q[0];
rz(-2.0149196) q[0];
sx q[0];
rz(-0.60420369) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2755865) q[2];
sx q[2];
rz(-1.3129932) q[2];
sx q[2];
rz(-0.38391963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4891543) q[1];
sx q[1];
rz(-1.592318) q[1];
sx q[1];
rz(-2.6179594) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0992364) q[3];
sx q[3];
rz(-2.5109249) q[3];
sx q[3];
rz(1.809465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.135123) q[2];
sx q[2];
rz(-2.3389356) q[2];
sx q[2];
rz(2.6378677) q[2];
rz(1.6040365) q[3];
sx q[3];
rz(-1.1732912) q[3];
sx q[3];
rz(-1.5776207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2876494) q[0];
sx q[0];
rz(-0.33677736) q[0];
sx q[0];
rz(-2.7119998) q[0];
rz(2.8751539) q[1];
sx q[1];
rz(-2.3315505) q[1];
sx q[1];
rz(2.0727167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42035139) q[0];
sx q[0];
rz(-1.9801894) q[0];
sx q[0];
rz(-2.4024525) q[0];
rz(-pi) q[1];
rz(-2.339459) q[2];
sx q[2];
rz(-1.7960227) q[2];
sx q[2];
rz(-0.98728131) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4539448) q[1];
sx q[1];
rz(-1.3140716) q[1];
sx q[1];
rz(0.4554847) q[1];
rz(-pi) q[2];
rz(-2.3924218) q[3];
sx q[3];
rz(-2.7862644) q[3];
sx q[3];
rz(-0.65566777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8916919) q[2];
sx q[2];
rz(-2.4155858) q[2];
sx q[2];
rz(2.967584) q[2];
rz(1.8788667) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(-1.5549829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1332909) q[0];
sx q[0];
rz(-1.0187848) q[0];
sx q[0];
rz(0.89023501) q[0];
rz(1.7091735) q[1];
sx q[1];
rz(-1.018367) q[1];
sx q[1];
rz(-0.11996809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9667119) q[0];
sx q[0];
rz(-0.57583664) q[0];
sx q[0];
rz(-2.8737646) q[0];
rz(-0.89996296) q[2];
sx q[2];
rz(-0.53027486) q[2];
sx q[2];
rz(-1.2880304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97291291) q[1];
sx q[1];
rz(-0.69818234) q[1];
sx q[1];
rz(-0.38777988) q[1];
rz(-pi) q[2];
rz(-2.8268485) q[3];
sx q[3];
rz(-2.5786244) q[3];
sx q[3];
rz(-1.986435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.031781901) q[2];
sx q[2];
rz(-1.6620771) q[2];
sx q[2];
rz(-0.39153448) q[2];
rz(1.3854965) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(0.63867205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93381298) q[0];
sx q[0];
rz(-2.5901828) q[0];
sx q[0];
rz(1.3054003) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.7326109) q[1];
sx q[1];
rz(-0.22444589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674057) q[0];
sx q[0];
rz(-1.5541557) q[0];
sx q[0];
rz(-1.5865302) q[0];
rz(-pi) q[1];
rz(-2.7573757) q[2];
sx q[2];
rz(-0.52500341) q[2];
sx q[2];
rz(2.8872761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.689381) q[1];
sx q[1];
rz(-1.2465917) q[1];
sx q[1];
rz(-0.28181847) q[1];
rz(-pi) q[2];
rz(-2.8185063) q[3];
sx q[3];
rz(-1.9226908) q[3];
sx q[3];
rz(-3.0754967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18181431) q[2];
sx q[2];
rz(-1.10428) q[2];
sx q[2];
rz(-2.1941198) q[2];
rz(-0.13253658) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(2.1555677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2333616) q[0];
sx q[0];
rz(-0.34603226) q[0];
sx q[0];
rz(-0.37286266) q[0];
rz(-2.8064959) q[1];
sx q[1];
rz(-1.2058039) q[1];
sx q[1];
rz(1.0710517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0308353) q[0];
sx q[0];
rz(-1.5016775) q[0];
sx q[0];
rz(1.5594307) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3577698) q[2];
sx q[2];
rz(-1.0292849) q[2];
sx q[2];
rz(1.0705186) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67642405) q[1];
sx q[1];
rz(-1.3260726) q[1];
sx q[1];
rz(-1.8108435) q[1];
x q[2];
rz(-1.2921807) q[3];
sx q[3];
rz(-2.8077586) q[3];
sx q[3];
rz(0.95164559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7143453) q[2];
sx q[2];
rz(-1.599396) q[2];
sx q[2];
rz(1.0672807) q[2];
rz(0.37211564) q[3];
sx q[3];
rz(-2.1477063) q[3];
sx q[3];
rz(2.9414162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2315955) q[0];
sx q[0];
rz(-0.59760004) q[0];
sx q[0];
rz(-1.7344612) q[0];
rz(1.0435957) q[1];
sx q[1];
rz(-2.2035972) q[1];
sx q[1];
rz(2.6503906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255913) q[0];
sx q[0];
rz(-1.859731) q[0];
sx q[0];
rz(-2.2641473) q[0];
rz(-pi) q[1];
rz(-1.4004854) q[2];
sx q[2];
rz(-1.404377) q[2];
sx q[2];
rz(-0.8569878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72043967) q[1];
sx q[1];
rz(-2.7758088) q[1];
sx q[1];
rz(0.37558742) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3902499) q[3];
sx q[3];
rz(-1.9420764) q[3];
sx q[3];
rz(0.89902395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4444943) q[2];
sx q[2];
rz(-1.631087) q[2];
sx q[2];
rz(1.0443045) q[2];
rz(0.77110428) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(2.7842298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412909) q[0];
sx q[0];
rz(-1.2933949) q[0];
sx q[0];
rz(1.6218761) q[0];
rz(-1.3007851) q[1];
sx q[1];
rz(-1.8404574) q[1];
sx q[1];
rz(0.69449743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.886777) q[0];
sx q[0];
rz(-0.22829994) q[0];
sx q[0];
rz(-1.6312902) q[0];
rz(-2.9143686) q[2];
sx q[2];
rz(-0.53582155) q[2];
sx q[2];
rz(0.69577755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3292103) q[1];
sx q[1];
rz(-0.33420104) q[1];
sx q[1];
rz(1.6045531) q[1];
rz(-pi) q[2];
rz(0.41037922) q[3];
sx q[3];
rz(-2.346092) q[3];
sx q[3];
rz(-2.8727227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0731611) q[2];
sx q[2];
rz(-1.3621) q[2];
sx q[2];
rz(2.6279602) q[2];
rz(0.3178151) q[3];
sx q[3];
rz(-1.603926) q[3];
sx q[3];
rz(-1.0191127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1401055) q[0];
sx q[0];
rz(-1.5167863) q[0];
sx q[0];
rz(2.232502) q[0];
rz(-1.9152676) q[1];
sx q[1];
rz(-1.5227804) q[1];
sx q[1];
rz(-2.7979122) q[1];
rz(-2.5334755) q[2];
sx q[2];
rz(-0.79610745) q[2];
sx q[2];
rz(2.983762) q[2];
rz(1.5671586) q[3];
sx q[3];
rz(-2.6156577) q[3];
sx q[3];
rz(-1.2729264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
