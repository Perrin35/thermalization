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
rz(-1.3353424) q[0];
sx q[0];
rz(3.5080533) q[0];
sx q[0];
rz(8.9656497) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(4.876457) q[1];
sx q[1];
rz(9.193037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9137751) q[0];
sx q[0];
rz(-1.639606) q[0];
sx q[0];
rz(-2.7976996) q[0];
x q[1];
rz(-2.5951573) q[2];
sx q[2];
rz(-1.6197512) q[2];
sx q[2];
rz(3.0374683) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79855483) q[1];
sx q[1];
rz(-0.75115582) q[1];
sx q[1];
rz(-2.5704774) q[1];
rz(-pi) q[2];
rz(3.0240293) q[3];
sx q[3];
rz(-2.8010578) q[3];
sx q[3];
rz(-1.3887608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0509402) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(1.862662) q[2];
rz(-0.88979641) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58562529) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(2.2157748) q[0];
rz(-2.0766808) q[1];
sx q[1];
rz(-1.5771259) q[1];
sx q[1];
rz(-0.085478641) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32498549) q[0];
sx q[0];
rz(-0.38982911) q[0];
sx q[0];
rz(3.1329324) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95909848) q[2];
sx q[2];
rz(-1.8597721) q[2];
sx q[2];
rz(-2.0069569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1053727) q[1];
sx q[1];
rz(-0.58507996) q[1];
sx q[1];
rz(1.9963422) q[1];
x q[2];
rz(-2.5898629) q[3];
sx q[3];
rz(-2.3750616) q[3];
sx q[3];
rz(-2.6006002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3257137) q[2];
sx q[2];
rz(-1.6625983) q[2];
sx q[2];
rz(0.386664) q[2];
rz(1.2169085) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(-2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3306408) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(2.6445342) q[0];
rz(-1.0423543) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(-0.29464468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2089068) q[0];
sx q[0];
rz(-2.8021376) q[0];
sx q[0];
rz(0.20821388) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63906901) q[2];
sx q[2];
rz(-1.4316171) q[2];
sx q[2];
rz(-2.7509909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5897419) q[1];
sx q[1];
rz(-2.183118) q[1];
sx q[1];
rz(0.53770868) q[1];
rz(-pi) q[2];
rz(-2.1717242) q[3];
sx q[3];
rz(-1.4517541) q[3];
sx q[3];
rz(0.62059072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7872539) q[2];
sx q[2];
rz(-0.86091176) q[2];
sx q[2];
rz(-1.5798689) q[2];
rz(2.0653557) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266066) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(2.9122747) q[0];
rz(-1.6353105) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(-0.062072676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631043) q[0];
sx q[0];
rz(-2.0651428) q[0];
sx q[0];
rz(-0.93017471) q[0];
rz(-0.27917316) q[2];
sx q[2];
rz(-1.9491157) q[2];
sx q[2];
rz(-2.407848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2885372) q[1];
sx q[1];
rz(-1.2934577) q[1];
sx q[1];
rz(-0.71373516) q[1];
rz(-2.9406406) q[3];
sx q[3];
rz(-1.4794599) q[3];
sx q[3];
rz(-0.82543711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.092992358) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(1.830706) q[2];
rz(0.038289573) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(2.766975) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467317) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(2.2485961) q[0];
rz(2.112174) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(0.095887862) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5928837) q[0];
sx q[0];
rz(-2.0796695) q[0];
sx q[0];
rz(0.48046099) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.65972) q[2];
sx q[2];
rz(-1.5006353) q[2];
sx q[2];
rz(1.905575) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8687415) q[1];
sx q[1];
rz(-0.68528803) q[1];
sx q[1];
rz(0.59602316) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60378243) q[3];
sx q[3];
rz(-0.60294916) q[3];
sx q[3];
rz(-2.3100738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9868077) q[2];
sx q[2];
rz(-1.3900577) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(2.7191539) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(-2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130037) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-2.2393176) q[1];
sx q[1];
rz(2.07043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18497047) q[0];
sx q[0];
rz(-1.8529467) q[0];
sx q[0];
rz(-0.50392242) q[0];
rz(2.5994634) q[2];
sx q[2];
rz(-1.8111992) q[2];
sx q[2];
rz(1.2803276) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6630604) q[1];
sx q[1];
rz(-2.2473865) q[1];
sx q[1];
rz(2.8512136) q[1];
rz(-pi) q[2];
rz(-0.46573205) q[3];
sx q[3];
rz(-0.35251401) q[3];
sx q[3];
rz(-2.5172212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57227197) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(0.039483698) q[2];
rz(-0.076400541) q[3];
sx q[3];
rz(-1.971784) q[3];
sx q[3];
rz(-0.45342818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385248) q[0];
sx q[0];
rz(-2.330133) q[0];
sx q[0];
rz(2.1912498) q[0];
rz(1.9901265) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(-1.8340402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81539736) q[0];
sx q[0];
rz(-1.895825) q[0];
sx q[0];
rz(2.9103878) q[0];
rz(-pi) q[1];
rz(1.6436395) q[2];
sx q[2];
rz(-1.9460287) q[2];
sx q[2];
rz(-0.90922395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89245172) q[1];
sx q[1];
rz(-1.080852) q[1];
sx q[1];
rz(2.9016414) q[1];
rz(-2.1767307) q[3];
sx q[3];
rz(-1.4601225) q[3];
sx q[3];
rz(-2.0697274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3211956) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(2.5879228) q[2];
rz(-0.6959483) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(-1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221136) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(2.8346862) q[0];
rz(1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(1.9006405) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0725482) q[0];
sx q[0];
rz(-2.6676151) q[0];
sx q[0];
rz(-1.9182253) q[0];
x q[1];
rz(-1.9996793) q[2];
sx q[2];
rz(-2.1253573) q[2];
sx q[2];
rz(0.51046023) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1859295) q[1];
sx q[1];
rz(-2.7765882) q[1];
sx q[1];
rz(-2.5899653) q[1];
rz(-0.60997529) q[3];
sx q[3];
rz(-1.2362567) q[3];
sx q[3];
rz(-1.9639249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3500195) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(-1.8512858) q[2];
rz(2.4604515) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(1.5017989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27555585) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(-1.1736897) q[0];
rz(-0.58700079) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(-0.079708286) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22686887) q[0];
sx q[0];
rz(-0.43926111) q[0];
sx q[0];
rz(2.125581) q[0];
rz(1.4585895) q[2];
sx q[2];
rz(-1.7908899) q[2];
sx q[2];
rz(2.8518563) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4533055) q[1];
sx q[1];
rz(-2.2548179) q[1];
sx q[1];
rz(2.3832393) q[1];
x q[2];
rz(0.2853197) q[3];
sx q[3];
rz(-0.82672182) q[3];
sx q[3];
rz(-0.43275012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6012663) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(1.5506844) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(-0.18868119) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48225668) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(1.2905066) q[0];
rz(3.105063) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(1.0443002) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14257763) q[0];
sx q[0];
rz(-1.8685088) q[0];
sx q[0];
rz(2.0172869) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9007334) q[2];
sx q[2];
rz(-1.3173027) q[2];
sx q[2];
rz(1.9869177) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83936939) q[1];
sx q[1];
rz(-1.7249134) q[1];
sx q[1];
rz(-0.83854143) q[1];
x q[2];
rz(2.812522) q[3];
sx q[3];
rz(-1.6794723) q[3];
sx q[3];
rz(1.4377126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9145987) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(-1.8168137) q[2];
rz(2.3850208) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-0.85048401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.666438) q[0];
sx q[0];
rz(-1.403724) q[0];
sx q[0];
rz(-2.4028461) q[0];
rz(-1.7186164) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(2.4074211) q[2];
sx q[2];
rz(-1.2470857) q[2];
sx q[2];
rz(-1.7294711) q[2];
rz(-2.041009) q[3];
sx q[3];
rz(-0.96405021) q[3];
sx q[3];
rz(-1.8586803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
