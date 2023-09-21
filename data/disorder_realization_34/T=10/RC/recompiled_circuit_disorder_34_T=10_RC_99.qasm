OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(0.08131942) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(4.4250017) q[1];
sx q[1];
rz(11.783574) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1845707) q[0];
sx q[0];
rz(-1.6324568) q[0];
sx q[0];
rz(1.8210568) q[0];
rz(1.9822257) q[2];
sx q[2];
rz(-1.4268268) q[2];
sx q[2];
rz(1.6978482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.767747) q[1];
sx q[1];
rz(-2.0031553) q[1];
sx q[1];
rz(-3.0004764) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8738535) q[3];
sx q[3];
rz(-0.32666884) q[3];
sx q[3];
rz(-2.3501345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.2540934) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-0.23981747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6759684) q[0];
sx q[0];
rz(-1.8166607) q[0];
sx q[0];
rz(-1.4093536) q[0];
rz(-pi) q[1];
rz(-2.1755978) q[2];
sx q[2];
rz(-2.7808393) q[2];
sx q[2];
rz(-0.15168562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0071734) q[1];
sx q[1];
rz(-1.6234142) q[1];
sx q[1];
rz(0.79018553) q[1];
x q[2];
rz(2.2858743) q[3];
sx q[3];
rz(-0.021589605) q[3];
sx q[3];
rz(1.74687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6372765) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(-1.3416393) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(-1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-0.66816107) q[0];
rz(-1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(1.0659165) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.802793) q[0];
sx q[0];
rz(-2.8899513) q[0];
sx q[0];
rz(-1.405707) q[0];
rz(-2.5862323) q[2];
sx q[2];
rz(-1.6931603) q[2];
sx q[2];
rz(-0.33081474) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6560116) q[1];
sx q[1];
rz(-2.1335568) q[1];
sx q[1];
rz(-0.33388168) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7245674) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9540017) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(0.70670635) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(-1.4867841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0006205) q[0];
sx q[0];
rz(-1.2396493) q[0];
sx q[0];
rz(-1.2739869) q[0];
rz(-pi) q[1];
rz(2.6629278) q[2];
sx q[2];
rz(-0.66648167) q[2];
sx q[2];
rz(-1.1043617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6390037) q[1];
sx q[1];
rz(-1.8559615) q[1];
sx q[1];
rz(-1.4147948) q[1];
x q[2];
rz(-1.5781457) q[3];
sx q[3];
rz(-0.74439936) q[3];
sx q[3];
rz(0.053754427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.460357) q[0];
rz(1.5530855) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(3.1255186) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094729312) q[0];
sx q[0];
rz(-1.0445347) q[0];
sx q[0];
rz(-0.88386436) q[0];
rz(0.0083382567) q[2];
sx q[2];
rz(-1.1401046) q[2];
sx q[2];
rz(-2.8991933) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4453455) q[1];
sx q[1];
rz(-1.0187341) q[1];
sx q[1];
rz(2.6389883) q[1];
x q[2];
rz(-0.36230476) q[3];
sx q[3];
rz(-1.5665896) q[3];
sx q[3];
rz(-1.6844974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(-2.0444929) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(0.90676701) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(1.9168568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46366102) q[0];
sx q[0];
rz(-1.7860798) q[0];
sx q[0];
rz(-0.93925516) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20829006) q[2];
sx q[2];
rz(-2.0049094) q[2];
sx q[2];
rz(-1.8722033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9310301) q[1];
sx q[1];
rz(-1.5157053) q[1];
sx q[1];
rz(0.964826) q[1];
rz(-1.46035) q[3];
sx q[3];
rz(-1.9697646) q[3];
sx q[3];
rz(-2.9810867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99888745) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9796824) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-2.5612223) q[0];
rz(-2.0866108) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(0.70077983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.472815) q[0];
sx q[0];
rz(-0.46572177) q[0];
sx q[0];
rz(-1.8229683) q[0];
rz(-pi) q[1];
rz(2.8220196) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(3.0802397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.339387) q[1];
sx q[1];
rz(-1.4805668) q[1];
sx q[1];
rz(2.1369364) q[1];
x q[2];
rz(3.0573781) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(-0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(-0.74470216) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3547524) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(1.7658866) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(-0.8917121) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211593) q[0];
sx q[0];
rz(-1.7283068) q[0];
sx q[0];
rz(-0.79173761) q[0];
x q[1];
rz(0.6363836) q[2];
sx q[2];
rz(-2.4578265) q[2];
sx q[2];
rz(0.46905876) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29854952) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(2.4119853) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1710839) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(1.2801998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(0.78424224) q[2];
rz(-2.6358321) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(-0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(1.3649712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46699539) q[0];
sx q[0];
rz(-2.31384) q[0];
sx q[0];
rz(-1.6270431) q[0];
rz(1.3806254) q[2];
sx q[2];
rz(-0.41751465) q[2];
sx q[2];
rz(-1.6585569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2862128) q[1];
sx q[1];
rz(-0.78821048) q[1];
sx q[1];
rz(-2.6469995) q[1];
x q[2];
rz(1.2737208) q[3];
sx q[3];
rz(-1.3798514) q[3];
sx q[3];
rz(1.2862658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(-1.2314679) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0868527) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(-2.1733984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4614842) q[0];
sx q[0];
rz(-2.619954) q[0];
sx q[0];
rz(-2.3011544) q[0];
rz(2.2553026) q[2];
sx q[2];
rz(-2.5286525) q[2];
sx q[2];
rz(-0.71000368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4149975) q[1];
sx q[1];
rz(-1.2508878) q[1];
sx q[1];
rz(1.725561) q[1];
rz(-0.00023437436) q[3];
sx q[3];
rz(-1.2486542) q[3];
sx q[3];
rz(-0.57516592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0633462) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(-2.3920849) q[2];
sx q[2];
rz(-1.8029677) q[2];
sx q[2];
rz(0.36063902) q[2];
rz(1.0352186) q[3];
sx q[3];
rz(-2.513701) q[3];
sx q[3];
rz(2.7660478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];