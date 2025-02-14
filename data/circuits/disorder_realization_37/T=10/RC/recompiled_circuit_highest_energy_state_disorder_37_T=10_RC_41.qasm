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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(-2.6824644) q[0];
rz(-0.65161172) q[1];
sx q[1];
rz(-1.7348644) q[1];
sx q[1];
rz(-2.9098517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2278175) q[0];
sx q[0];
rz(-1.5019866) q[0];
sx q[0];
rz(2.7976996) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5135147) q[2];
sx q[2];
rz(-2.1165032) q[2];
sx q[2];
rz(1.4964262) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3430378) q[1];
sx q[1];
rz(-2.3904368) q[1];
sx q[1];
rz(0.57111528) q[1];
rz(-1.5292589) q[3];
sx q[3];
rz(-1.2327063) q[3];
sx q[3];
rz(-1.6281782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0509402) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(-1.862662) q[2];
rz(0.88979641) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(-0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58562529) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(-0.92581785) q[0];
rz(-1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(3.056114) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8259698) q[0];
sx q[0];
rz(-1.1809826) q[0];
sx q[0];
rz(1.5743544) q[0];
rz(0.95909848) q[2];
sx q[2];
rz(-1.2818205) q[2];
sx q[2];
rz(-2.0069569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6798415) q[1];
sx q[1];
rz(-2.0978758) q[1];
sx q[1];
rz(-0.26694571) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5898629) q[3];
sx q[3];
rz(-0.76653102) q[3];
sx q[3];
rz(2.6006002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3257137) q[2];
sx q[2];
rz(-1.6625983) q[2];
sx q[2];
rz(2.7549287) q[2];
rz(1.9246842) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(-2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3306408) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(0.49705848) q[0];
rz(-1.0423543) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(-2.846948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9884856) q[0];
sx q[0];
rz(-1.238958) q[0];
sx q[0];
rz(-1.4979304) q[0];
rz(-pi) q[1];
rz(-2.5025236) q[2];
sx q[2];
rz(-1.7099755) q[2];
sx q[2];
rz(-0.39060171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5897419) q[1];
sx q[1];
rz(-2.183118) q[1];
sx q[1];
rz(-0.53770868) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96986846) q[3];
sx q[3];
rz(-1.6898385) q[3];
sx q[3];
rz(2.5210019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3543388) q[2];
sx q[2];
rz(-0.86091176) q[2];
sx q[2];
rz(-1.5617237) q[2];
rz(-2.0653557) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(-2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.314986) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(2.9122747) q[0];
rz(-1.6353105) q[1];
sx q[1];
rz(-1.458026) q[1];
sx q[1];
rz(-3.07952) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45907606) q[0];
sx q[0];
rz(-2.354265) q[0];
sx q[0];
rz(-2.3045901) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1770085) q[2];
sx q[2];
rz(-2.6754489) q[2];
sx q[2];
rz(0.073748253) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.049438795) q[1];
sx q[1];
rz(-2.2518932) q[1];
sx q[1];
rz(-1.9309631) q[1];
rz(-pi) q[2];
rz(0.43020474) q[3];
sx q[3];
rz(-2.9211126) q[3];
sx q[3];
rz(1.9752432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0486003) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(-1.3108866) q[2];
rz(3.1033031) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(-2.766975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49486092) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(-2.2485961) q[0];
rz(-2.112174) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(-0.095887862) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9150118) q[0];
sx q[0];
rz(-1.1553197) q[0];
sx q[0];
rz(-1.0092495) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0711543) q[2];
sx q[2];
rz(-1.6595006) q[2];
sx q[2];
rz(2.8005637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27285114) q[1];
sx q[1];
rz(-0.68528803) q[1];
sx q[1];
rz(0.59602316) q[1];
x q[2];
rz(-0.60378243) q[3];
sx q[3];
rz(-2.5386435) q[3];
sx q[3];
rz(2.3100738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9868077) q[2];
sx q[2];
rz(-1.3900577) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(2.6384242) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4130037) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(3.0244306) q[0];
rz(1.4940184) q[1];
sx q[1];
rz(-2.2393176) q[1];
sx q[1];
rz(1.0711627) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9566222) q[0];
sx q[0];
rz(-1.288646) q[0];
sx q[0];
rz(0.50392242) q[0];
rz(0.54212924) q[2];
sx q[2];
rz(-1.8111992) q[2];
sx q[2];
rz(1.8612651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92381645) q[1];
sx q[1];
rz(-2.4144396) q[1];
sx q[1];
rz(-1.2283065) q[1];
rz(0.31757327) q[3];
sx q[3];
rz(-1.7264719) q[3];
sx q[3];
rz(-0.50567108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57227197) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(-0.039483698) q[2];
rz(3.0651921) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(0.45342818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-2.1912498) q[0];
rz(-1.1514661) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(1.8340402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83043419) q[0];
sx q[0];
rz(-1.7896928) q[0];
sx q[0];
rz(1.9040742) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7654539) q[2];
sx q[2];
rz(-1.6385632) q[2];
sx q[2];
rz(2.4532831) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7286018) q[1];
sx q[1];
rz(-2.6003693) q[1];
sx q[1];
rz(1.9899998) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96486196) q[3];
sx q[3];
rz(-1.4601225) q[3];
sx q[3];
rz(-2.0697274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82039708) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(2.5879228) q[2];
rz(2.4456444) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(1.6863916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11947908) q[0];
sx q[0];
rz(-1.6895634) q[0];
sx q[0];
rz(2.8346862) q[0];
rz(1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(-1.2409522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3281455) q[0];
sx q[0];
rz(-1.4147583) q[0];
sx q[0];
rz(-1.1213844) q[0];
rz(-0.5979171) q[2];
sx q[2];
rz(-1.9321402) q[2];
sx q[2];
rz(1.2966228) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7684264) q[1];
sx q[1];
rz(-1.2618999) q[1];
sx q[1];
rz(1.7684446) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5316174) q[3];
sx q[3];
rz(-1.2362567) q[3];
sx q[3];
rz(-1.1776678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79157311) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(-1.2903068) q[2];
rz(-2.4604515) q[3];
sx q[3];
rz(-0.9001503) q[3];
sx q[3];
rz(-1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27555585) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(-1.9679029) q[0];
rz(-2.5545919) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(3.0618844) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3088206) q[0];
sx q[0];
rz(-1.7967294) q[0];
sx q[0];
rz(-1.1908047) q[0];
rz(1.4585895) q[2];
sx q[2];
rz(-1.7908899) q[2];
sx q[2];
rz(2.8518563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4707798) q[1];
sx q[1];
rz(-0.97320405) q[1];
sx q[1];
rz(0.87009354) q[1];
rz(-pi) q[2];
rz(0.80612889) q[3];
sx q[3];
rz(-1.7793831) q[3];
sx q[3];
rz(-0.94193469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54032636) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(1.5909083) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(-2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659336) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(1.2905066) q[0];
rz(0.036529649) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(2.0972924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8528906) q[0];
sx q[0];
rz(-1.1452617) q[0];
sx q[0];
rz(-0.3279) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89628156) q[2];
sx q[2];
rz(-2.7283629) q[2];
sx q[2];
rz(-0.21597029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56230169) q[1];
sx q[1];
rz(-2.3962483) q[1];
sx q[1];
rz(1.3424804) q[1];
rz(-1.6855816) q[3];
sx q[3];
rz(-1.8978531) q[3];
sx q[3];
rz(-2.9714874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(1.3247789) q[2];
rz(-2.3850208) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-2.2911086) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751547) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(1.7186164) q[1];
sx q[1];
rz(-1.9444793) q[1];
sx q[1];
rz(0.3107298) q[1];
rz(-2.4074211) q[2];
sx q[2];
rz(-1.894507) q[2];
sx q[2];
rz(1.4121216) q[2];
rz(-0.66154578) q[3];
sx q[3];
rz(-1.9521803) q[3];
sx q[3];
rz(3.1357756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
