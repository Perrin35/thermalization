OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(2.0318883) q[0];
sx q[0];
rz(11.561031) q[0];
rz(0.73468626) q[1];
sx q[1];
rz(3.4975657) q[1];
sx q[1];
rz(11.674292) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1901613) q[0];
sx q[0];
rz(-2.3888458) q[0];
sx q[0];
rz(-2.997082) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1744035) q[2];
sx q[2];
rz(-1.4481067) q[2];
sx q[2];
rz(-1.762378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0301117) q[1];
sx q[1];
rz(-1.6684384) q[1];
sx q[1];
rz(0.48508118) q[1];
rz(-1.0451116) q[3];
sx q[3];
rz(-1.0307923) q[3];
sx q[3];
rz(-1.2837376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72508183) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(1.7007281) q[2];
rz(2.1848988) q[3];
sx q[3];
rz(-1.1172833) q[3];
sx q[3];
rz(-2.8982437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6956534) q[0];
sx q[0];
rz(-2.74701) q[0];
sx q[0];
rz(1.12895) q[0];
rz(2.8961862) q[1];
sx q[1];
rz(-1.3134198) q[1];
sx q[1];
rz(0.66171563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11479325) q[0];
sx q[0];
rz(-0.68959348) q[0];
sx q[0];
rz(-0.45701509) q[0];
x q[1];
rz(1.7956338) q[2];
sx q[2];
rz(-1.1974466) q[2];
sx q[2];
rz(0.38885798) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0894818) q[1];
sx q[1];
rz(-2.8683337) q[1];
sx q[1];
rz(-0.11788003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8850967) q[3];
sx q[3];
rz(-0.30499015) q[3];
sx q[3];
rz(1.631032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5976065) q[2];
sx q[2];
rz(-0.49763766) q[2];
sx q[2];
rz(2.0987161) q[2];
rz(1.4526224) q[3];
sx q[3];
rz(-2.2904604) q[3];
sx q[3];
rz(-1.1037306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76688981) q[0];
sx q[0];
rz(-2.4536528) q[0];
sx q[0];
rz(1.8324628) q[0];
rz(0.4492999) q[1];
sx q[1];
rz(-2.4520912) q[1];
sx q[1];
rz(1.4264533) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5912522) q[0];
sx q[0];
rz(-0.36446291) q[0];
sx q[0];
rz(2.3191207) q[0];
x q[1];
rz(-1.1482936) q[2];
sx q[2];
rz(-1.587217) q[2];
sx q[2];
rz(-0.33994477) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5876679) q[1];
sx q[1];
rz(-1.4933506) q[1];
sx q[1];
rz(-0.70550349) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3247812) q[3];
sx q[3];
rz(-1.4895456) q[3];
sx q[3];
rz(1.2556835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7766777) q[2];
sx q[2];
rz(-1.1740351) q[2];
sx q[2];
rz(3.0752227) q[2];
rz(0.55473173) q[3];
sx q[3];
rz(-1.4812255) q[3];
sx q[3];
rz(-1.6248645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1944815) q[0];
sx q[0];
rz(-0.28488657) q[0];
sx q[0];
rz(2.3696005) q[0];
rz(-2.4783065) q[1];
sx q[1];
rz(-0.74378496) q[1];
sx q[1];
rz(3.0139121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049162) q[0];
sx q[0];
rz(-0.66670115) q[0];
sx q[0];
rz(-0.0078860869) q[0];
rz(-pi) q[1];
rz(-0.56945412) q[2];
sx q[2];
rz(-0.63210154) q[2];
sx q[2];
rz(-2.5394627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3290594) q[1];
sx q[1];
rz(-0.61186463) q[1];
sx q[1];
rz(-2.3298323) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8074266) q[3];
sx q[3];
rz(-0.45229707) q[3];
sx q[3];
rz(2.3826007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6105303) q[2];
sx q[2];
rz(-2.1805111) q[2];
sx q[2];
rz(0.5985716) q[2];
rz(-2.5564204) q[3];
sx q[3];
rz(-1.3734615) q[3];
sx q[3];
rz(-3.0857871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7912306) q[0];
sx q[0];
rz(-1.5086011) q[0];
sx q[0];
rz(-2.1685261) q[0];
rz(-2.9452501) q[1];
sx q[1];
rz(-2.0025496) q[1];
sx q[1];
rz(-1.4686718) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40590433) q[0];
sx q[0];
rz(-2.2668512) q[0];
sx q[0];
rz(0.52175867) q[0];
rz(-1.0650915) q[2];
sx q[2];
rz(-1.7943766) q[2];
sx q[2];
rz(0.93709842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10188516) q[1];
sx q[1];
rz(-0.3488003) q[1];
sx q[1];
rz(-1.0451911) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0861868) q[3];
sx q[3];
rz(-2.1089777) q[3];
sx q[3];
rz(-2.2937867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68957442) q[2];
sx q[2];
rz(-1.8147899) q[2];
sx q[2];
rz(-0.16239521) q[2];
rz(1.1299805) q[3];
sx q[3];
rz(-2.2584848) q[3];
sx q[3];
rz(1.1348772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75055403) q[0];
sx q[0];
rz(-2.6226608) q[0];
sx q[0];
rz(-3.0911875) q[0];
rz(-2.9734036) q[1];
sx q[1];
rz(-1.5805809) q[1];
sx q[1];
rz(2.2680297) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024603) q[0];
sx q[0];
rz(-0.96193571) q[0];
sx q[0];
rz(1.2723421) q[0];
rz(-1.7123772) q[2];
sx q[2];
rz(-1.537286) q[2];
sx q[2];
rz(0.17591116) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2312886) q[1];
sx q[1];
rz(-1.043348) q[1];
sx q[1];
rz(-2.0922679) q[1];
rz(-pi) q[2];
rz(-0.29804067) q[3];
sx q[3];
rz(-1.3750769) q[3];
sx q[3];
rz(-2.6258371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1217338) q[2];
sx q[2];
rz(-2.9254318) q[2];
sx q[2];
rz(-0.38062322) q[2];
rz(0.56626433) q[3];
sx q[3];
rz(-1.7411313) q[3];
sx q[3];
rz(-0.80999058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.64067632) q[0];
sx q[0];
rz(-0.2114507) q[0];
sx q[0];
rz(0.40263116) q[0];
rz(-1.3817878) q[1];
sx q[1];
rz(-0.80843061) q[1];
sx q[1];
rz(-0.056338739) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9026109) q[0];
sx q[0];
rz(-1.9704116) q[0];
sx q[0];
rz(0.35524551) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4903421) q[2];
sx q[2];
rz(-2.6695996) q[2];
sx q[2];
rz(2.9126963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.64526865) q[1];
sx q[1];
rz(-2.5773415) q[1];
sx q[1];
rz(1.5936046) q[1];
x q[2];
rz(-1.0545441) q[3];
sx q[3];
rz(-2.1813381) q[3];
sx q[3];
rz(0.96585376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9075883) q[2];
sx q[2];
rz(-1.3434709) q[2];
sx q[2];
rz(2.6410356) q[2];
rz(0.060700011) q[3];
sx q[3];
rz(-1.7874291) q[3];
sx q[3];
rz(-0.48262706) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57182264) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(-1.7806336) q[0];
rz(0.743615) q[1];
sx q[1];
rz(-1.7230325) q[1];
sx q[1];
rz(-0.13882151) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25538975) q[0];
sx q[0];
rz(-2.3950044) q[0];
sx q[0];
rz(-2.3476178) q[0];
rz(0.93075625) q[2];
sx q[2];
rz(-2.4413707) q[2];
sx q[2];
rz(2.8559358) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.946859) q[1];
sx q[1];
rz(-1.3689965) q[1];
sx q[1];
rz(1.9029593) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0431792) q[3];
sx q[3];
rz(-2.7076027) q[3];
sx q[3];
rz(-1.7349617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4678141) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(-0.71411258) q[2];
rz(2.1821187) q[3];
sx q[3];
rz(-1.4941314) q[3];
sx q[3];
rz(0.97904557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7239083) q[0];
sx q[0];
rz(-2.4651616) q[0];
sx q[0];
rz(-3.0133001) q[0];
rz(2.7543606) q[1];
sx q[1];
rz(-0.59806824) q[1];
sx q[1];
rz(1.8181575) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5538226) q[0];
sx q[0];
rz(-2.9284796) q[0];
sx q[0];
rz(-2.0487259) q[0];
rz(0.88487423) q[2];
sx q[2];
rz(-2.0506244) q[2];
sx q[2];
rz(-0.72140933) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7625757) q[1];
sx q[1];
rz(-1.5847995) q[1];
sx q[1];
rz(-1.714731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72800379) q[3];
sx q[3];
rz(-1.4211072) q[3];
sx q[3];
rz(2.585175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4667751) q[2];
sx q[2];
rz(-1.520949) q[2];
sx q[2];
rz(0.20469323) q[2];
rz(1.2734867) q[3];
sx q[3];
rz(-2.4114362) q[3];
sx q[3];
rz(1.5654806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8650763) q[0];
sx q[0];
rz(-0.62224329) q[0];
sx q[0];
rz(2.9283071) q[0];
rz(0.22008303) q[1];
sx q[1];
rz(-1.8890231) q[1];
sx q[1];
rz(-1.782104) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3228139) q[0];
sx q[0];
rz(-1.5548163) q[0];
sx q[0];
rz(-3.131991) q[0];
rz(2.67138) q[2];
sx q[2];
rz(-0.5222975) q[2];
sx q[2];
rz(0.14005113) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46392808) q[1];
sx q[1];
rz(-1.0654133) q[1];
sx q[1];
rz(2.7457854) q[1];
rz(-pi) q[2];
rz(2.241019) q[3];
sx q[3];
rz(-1.9346721) q[3];
sx q[3];
rz(-2.8577141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6243817) q[2];
sx q[2];
rz(-1.5959847) q[2];
sx q[2];
rz(0.33779302) q[2];
rz(1.4126011) q[3];
sx q[3];
rz(-1.9905041) q[3];
sx q[3];
rz(-0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4008041) q[0];
sx q[0];
rz(-0.82397006) q[0];
sx q[0];
rz(-1.9116221) q[0];
rz(-0.38372718) q[1];
sx q[1];
rz(-1.6576672) q[1];
sx q[1];
rz(-2.3794649) q[1];
rz(-0.40390759) q[2];
sx q[2];
rz(-1.8941034) q[2];
sx q[2];
rz(-0.20878172) q[2];
rz(0.2445515) q[3];
sx q[3];
rz(-2.0691732) q[3];
sx q[3];
rz(0.3680784) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
