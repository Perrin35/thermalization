OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95028967) q[0];
sx q[0];
rz(6.0130881) q[0];
sx q[0];
rz(10.313378) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(-1.5421901) q[1];
sx q[1];
rz(1.3808274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61533538) q[0];
sx q[0];
rz(-1.3852296) q[0];
sx q[0];
rz(-2.944817) q[0];
rz(0.51989748) q[2];
sx q[2];
rz(-2.2714104) q[2];
sx q[2];
rz(-1.1288647) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12871582) q[1];
sx q[1];
rz(-2.8059462) q[1];
sx q[1];
rz(-0.0291834) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9350151) q[3];
sx q[3];
rz(-1.7813588) q[3];
sx q[3];
rz(-0.59277804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8218653) q[2];
sx q[2];
rz(-1.8165908) q[2];
sx q[2];
rz(-0.74903178) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(-0.099893959) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(-1.9248167) q[0];
rz(1.0617537) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(-0.42713508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9574643) q[0];
sx q[0];
rz(-1.2250568) q[0];
sx q[0];
rz(0.60681245) q[0];
rz(-pi) q[1];
rz(0.92556503) q[2];
sx q[2];
rz(-1.3750018) q[2];
sx q[2];
rz(-2.2044971) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9717025) q[1];
sx q[1];
rz(-2.8198543) q[1];
sx q[1];
rz(1.1119026) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8187567) q[3];
sx q[3];
rz(-2.9313847) q[3];
sx q[3];
rz(1.5402286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(1.742935) q[2];
rz(0.22377293) q[3];
sx q[3];
rz(-2.2327435) q[3];
sx q[3];
rz(1.5980501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.021521213) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(-0.045510005) q[0];
rz(-1.1687763) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(0.57560903) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110014) q[0];
sx q[0];
rz(-1.5819893) q[0];
sx q[0];
rz(-2.5024947) q[0];
rz(-pi) q[1];
rz(-1.5928629) q[2];
sx q[2];
rz(-2.4808072) q[2];
sx q[2];
rz(1.929785) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2896881) q[1];
sx q[1];
rz(-0.59483268) q[1];
sx q[1];
rz(-1.1762876) q[1];
rz(-pi) q[2];
rz(0.91797249) q[3];
sx q[3];
rz(-1.550972) q[3];
sx q[3];
rz(1.6719847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1316954) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(0.95727813) q[2];
rz(-0.57473985) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(1.3055698) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15605536) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(1.8804469) q[0];
rz(0.082322923) q[1];
sx q[1];
rz(-2.0539093) q[1];
sx q[1];
rz(-1.3538768) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.026431) q[0];
sx q[0];
rz(-1.3360093) q[0];
sx q[0];
rz(0.26682667) q[0];
rz(-pi) q[1];
rz(1.3867845) q[2];
sx q[2];
rz(-1.0920326) q[2];
sx q[2];
rz(1.1637853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17300837) q[1];
sx q[1];
rz(-2.6436596) q[1];
sx q[1];
rz(-0.11346101) q[1];
rz(-pi) q[2];
rz(2.2500751) q[3];
sx q[3];
rz(-1.1783021) q[3];
sx q[3];
rz(0.52549997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7202619) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(1.8850231) q[2];
rz(-2.4300857) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8307777) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(-0.34314439) q[0];
rz(-3.0798196) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(-1.7247346) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8648711) q[0];
sx q[0];
rz(-1.7407932) q[0];
sx q[0];
rz(-1.8525339) q[0];
rz(-pi) q[1];
rz(0.38133905) q[2];
sx q[2];
rz(-0.73753192) q[2];
sx q[2];
rz(-0.87070891) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7570863) q[1];
sx q[1];
rz(-2.050274) q[1];
sx q[1];
rz(0.15285413) q[1];
rz(-pi) q[2];
rz(0.23483488) q[3];
sx q[3];
rz(-1.1560625) q[3];
sx q[3];
rz(2.107634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5477649) q[2];
sx q[2];
rz(-1.4676899) q[2];
sx q[2];
rz(1.0859547) q[2];
rz(0.91935277) q[3];
sx q[3];
rz(-1.7707526) q[3];
sx q[3];
rz(-0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.74349657) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(-2.3486775) q[0];
rz(-0.74053699) q[1];
sx q[1];
rz(-2.1426327) q[1];
sx q[1];
rz(0.83121306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8342469) q[0];
sx q[0];
rz(-1.9458658) q[0];
sx q[0];
rz(2.9177279) q[0];
rz(-0.0074499091) q[2];
sx q[2];
rz(-1.8478571) q[2];
sx q[2];
rz(-1.7817093) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97396353) q[1];
sx q[1];
rz(-2.2669753) q[1];
sx q[1];
rz(-2.7401398) q[1];
x q[2];
rz(1.6115723) q[3];
sx q[3];
rz(-2.6705461) q[3];
sx q[3];
rz(-0.99179964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.071216019) q[2];
sx q[2];
rz(-1.8698317) q[2];
sx q[2];
rz(-2.0224723) q[2];
rz(-1.2707155) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(-1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374461) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(1.5995837) q[0];
rz(1.0150602) q[1];
sx q[1];
rz(-1.5856182) q[1];
sx q[1];
rz(0.17280811) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5046523) q[0];
sx q[0];
rz(-2.0604134) q[0];
sx q[0];
rz(2.2235653) q[0];
rz(-pi) q[1];
rz(-2.6989231) q[2];
sx q[2];
rz(-0.73390642) q[2];
sx q[2];
rz(-1.1709605) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0902675) q[1];
sx q[1];
rz(-2.6410612) q[1];
sx q[1];
rz(1.0330908) q[1];
rz(-pi) q[2];
rz(-2.5739848) q[3];
sx q[3];
rz(-2.2877573) q[3];
sx q[3];
rz(-2.6102118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66199866) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(0.97770989) q[2];
rz(-2.5665723) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(1.9112816) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8136895) q[0];
sx q[0];
rz(-2.8459025) q[0];
sx q[0];
rz(-0.015722474) q[0];
rz(-1.9533336) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(-1.1994919) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6576783) q[0];
sx q[0];
rz(-1.8059314) q[0];
sx q[0];
rz(1.1170618) q[0];
rz(-pi) q[1];
x q[1];
rz(1.079915) q[2];
sx q[2];
rz(-1.617031) q[2];
sx q[2];
rz(0.85109988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53912341) q[1];
sx q[1];
rz(-2.0401461) q[1];
sx q[1];
rz(-0.55901171) q[1];
rz(-pi) q[2];
rz(-1.1114798) q[3];
sx q[3];
rz(-2.5190341) q[3];
sx q[3];
rz(0.42296577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1477995) q[2];
sx q[2];
rz(-1.2387929) q[2];
sx q[2];
rz(-1.7564868) q[2];
rz(0.6238474) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(2.0417716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1856336) q[0];
sx q[0];
rz(-0.82056844) q[0];
sx q[0];
rz(2.4483335) q[0];
rz(0.26501003) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(-1.5230491) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6054305) q[0];
sx q[0];
rz(-2.3443065) q[0];
sx q[0];
rz(1.5553586) q[0];
rz(0.9829282) q[2];
sx q[2];
rz(-2.1097906) q[2];
sx q[2];
rz(-1.3608152) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56227389) q[1];
sx q[1];
rz(-2.8011311) q[1];
sx q[1];
rz(3.0698983) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.070799) q[3];
sx q[3];
rz(-2.1507743) q[3];
sx q[3];
rz(1.1759315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8388464) q[2];
sx q[2];
rz(-1.5134209) q[2];
sx q[2];
rz(1.6839074) q[2];
rz(2.1863106) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38462287) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(-0.48267522) q[0];
rz(-0.92974281) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(0.39628705) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.020655) q[0];
sx q[0];
rz(-1.2115941) q[0];
sx q[0];
rz(0.44393702) q[0];
rz(-0.70102923) q[2];
sx q[2];
rz(-0.23244474) q[2];
sx q[2];
rz(1.6384517) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5037646) q[1];
sx q[1];
rz(-0.85496584) q[1];
sx q[1];
rz(2.9010335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40376264) q[3];
sx q[3];
rz(-1.9039394) q[3];
sx q[3];
rz(-3.1333095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79260176) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(2.8144042) q[2];
rz(2.3606825) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1031716) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(0.36956638) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(1.6160028) q[2];
sx q[2];
rz(-0.87799413) q[2];
sx q[2];
rz(-2.9064657) q[2];
rz(1.0450324) q[3];
sx q[3];
rz(-0.67787328) q[3];
sx q[3];
rz(-1.9280435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
