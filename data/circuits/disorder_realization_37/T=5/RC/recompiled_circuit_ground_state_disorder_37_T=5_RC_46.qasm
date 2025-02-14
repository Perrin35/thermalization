OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6310298) q[0];
sx q[0];
rz(3.6462311) q[0];
sx q[0];
rz(12.991821) q[0];
rz(-2.5692441) q[1];
sx q[1];
rz(-1.0971789) q[1];
sx q[1];
rz(-1.2001349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6365189) q[0];
sx q[0];
rz(-0.25101343) q[0];
sx q[0];
rz(-0.039029718) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91325735) q[2];
sx q[2];
rz(-2.7139671) q[2];
sx q[2];
rz(1.0853801) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26781005) q[1];
sx q[1];
rz(-1.7846196) q[1];
sx q[1];
rz(0.38114433) q[1];
rz(-pi) q[2];
rz(0.17961924) q[3];
sx q[3];
rz(-0.67905513) q[3];
sx q[3];
rz(1.2509026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9878865) q[2];
sx q[2];
rz(-1.989863) q[2];
sx q[2];
rz(-2.494334) q[2];
rz(0.17793947) q[3];
sx q[3];
rz(-1.2578332) q[3];
sx q[3];
rz(1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4030289) q[0];
sx q[0];
rz(-3.0572427) q[0];
sx q[0];
rz(-2.6179598) q[0];
rz(1.2298443) q[1];
sx q[1];
rz(-0.74408999) q[1];
sx q[1];
rz(-0.24807608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46892525) q[0];
sx q[0];
rz(-1.4644196) q[0];
sx q[0];
rz(2.399978) q[0];
rz(-pi) q[1];
rz(2.7047728) q[2];
sx q[2];
rz(-1.9736276) q[2];
sx q[2];
rz(-0.4536597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7404799) q[1];
sx q[1];
rz(-1.7164299) q[1];
sx q[1];
rz(3.1023953) q[1];
x q[2];
rz(2.5026191) q[3];
sx q[3];
rz(-0.77497122) q[3];
sx q[3];
rz(-2.4955622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6123885) q[2];
sx q[2];
rz(-0.91823053) q[2];
sx q[2];
rz(0.53981346) q[2];
rz(-2.9811033) q[3];
sx q[3];
rz(-1.5925946) q[3];
sx q[3];
rz(-0.81015712) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983343) q[0];
sx q[0];
rz(-2.46471) q[0];
sx q[0];
rz(-0.13370378) q[0];
rz(1.5191822) q[1];
sx q[1];
rz(-1.8459873) q[1];
sx q[1];
rz(-2.349283) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2871029) q[0];
sx q[0];
rz(-2.7407795) q[0];
sx q[0];
rz(1.0846179) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73826628) q[2];
sx q[2];
rz(-2.9472137) q[2];
sx q[2];
rz(-0.76972085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49323248) q[1];
sx q[1];
rz(-0.91714749) q[1];
sx q[1];
rz(-2.3834989) q[1];
x q[2];
rz(-1.5709223) q[3];
sx q[3];
rz(-1.2706666) q[3];
sx q[3];
rz(-2.2430994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2866659) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(2.3112042) q[2];
rz(1.7928127) q[3];
sx q[3];
rz(-2.1068137) q[3];
sx q[3];
rz(-0.59876284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778359) q[0];
sx q[0];
rz(-0.97665518) q[0];
sx q[0];
rz(-0.42309716) q[0];
rz(-1.1571723) q[1];
sx q[1];
rz(-1.4856228) q[1];
sx q[1];
rz(-0.62209904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4870479) q[0];
sx q[0];
rz(-1.9584987) q[0];
sx q[0];
rz(2.2564006) q[0];
rz(-pi) q[1];
rz(-1.9900722) q[2];
sx q[2];
rz(-1.7618693) q[2];
sx q[2];
rz(-2.9951468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8922975) q[1];
sx q[1];
rz(-1.163332) q[1];
sx q[1];
rz(-3.1235126) q[1];
rz(-pi) q[2];
rz(-1.5399718) q[3];
sx q[3];
rz(-0.72576952) q[3];
sx q[3];
rz(-2.9963065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28079924) q[2];
sx q[2];
rz(-1.4192899) q[2];
sx q[2];
rz(0.61817509) q[2];
rz(-1.6587967) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16267714) q[0];
sx q[0];
rz(-1.3293043) q[0];
sx q[0];
rz(0.83258122) q[0];
rz(2.816448) q[1];
sx q[1];
rz(-0.99601662) q[1];
sx q[1];
rz(-3.0772193) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20124931) q[0];
sx q[0];
rz(-1.8953249) q[0];
sx q[0];
rz(-1.6973116) q[0];
x q[1];
rz(0.60466296) q[2];
sx q[2];
rz(-0.69398601) q[2];
sx q[2];
rz(-2.1460905) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2657703) q[1];
sx q[1];
rz(-2.0811305) q[1];
sx q[1];
rz(0.93632621) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1983002) q[3];
sx q[3];
rz(-1.2767013) q[3];
sx q[3];
rz(2.7315745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9084106) q[2];
sx q[2];
rz(-2.5711214) q[2];
sx q[2];
rz(-1.849966) q[2];
rz(2.6455247) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(1.9991649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9555776) q[0];
sx q[0];
rz(-0.29009524) q[0];
sx q[0];
rz(0.41906038) q[0];
rz(-2.8298607) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(-0.14351235) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18565047) q[0];
sx q[0];
rz(-0.5910735) q[0];
sx q[0];
rz(-0.67703621) q[0];
rz(-1.8060051) q[2];
sx q[2];
rz(-1.303788) q[2];
sx q[2];
rz(0.63797229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.011605) q[1];
sx q[1];
rz(-1.2530348) q[1];
sx q[1];
rz(0.77316534) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49893283) q[3];
sx q[3];
rz(-0.29424516) q[3];
sx q[3];
rz(-0.1732451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.781337) q[2];
sx q[2];
rz(-2.2402253) q[2];
sx q[2];
rz(2.0085013) q[2];
rz(-1.1059149) q[3];
sx q[3];
rz(-1.7224885) q[3];
sx q[3];
rz(-2.6667986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7853506) q[0];
sx q[0];
rz(-2.3478735) q[0];
sx q[0];
rz(-1.445795) q[0];
rz(0.71714199) q[1];
sx q[1];
rz(-1.278806) q[1];
sx q[1];
rz(0.76146567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4091771) q[0];
sx q[0];
rz(-1.4237561) q[0];
sx q[0];
rz(0.63296403) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0324208) q[2];
sx q[2];
rz(-0.79029492) q[2];
sx q[2];
rz(0.99144713) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.059342) q[1];
sx q[1];
rz(-2.5656156) q[1];
sx q[1];
rz(-1.4800998) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7588571) q[3];
sx q[3];
rz(-1.4131964) q[3];
sx q[3];
rz(0.068305123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36934272) q[2];
sx q[2];
rz(-2.6476634) q[2];
sx q[2];
rz(-0.93144766) q[2];
rz(2.5233968) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(-2.1759694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44929993) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(-1.8071254) q[0];
rz(-0.5131228) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(-1.9122874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72237678) q[0];
sx q[0];
rz(-0.88787365) q[0];
sx q[0];
rz(2.1238668) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0487829) q[2];
sx q[2];
rz(-1.3121737) q[2];
sx q[2];
rz(1.5802204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7569983) q[1];
sx q[1];
rz(-1.726686) q[1];
sx q[1];
rz(-1.5499359) q[1];
rz(-pi) q[2];
rz(-3.119) q[3];
sx q[3];
rz(-0.66231662) q[3];
sx q[3];
rz(-2.3316771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3076155) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(-2.9311467) q[2];
rz(3.0757507) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(0.79894799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54302067) q[0];
sx q[0];
rz(-2.4871171) q[0];
sx q[0];
rz(1.2731113) q[0];
rz(-1.2741362) q[1];
sx q[1];
rz(-2.4576063) q[1];
sx q[1];
rz(0.88409105) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40128517) q[0];
sx q[0];
rz(-1.570556) q[0];
sx q[0];
rz(-0.9539286) q[0];
x q[1];
rz(2.6173148) q[2];
sx q[2];
rz(-1.345621) q[2];
sx q[2];
rz(2.2543627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6581304) q[1];
sx q[1];
rz(-2.6331055) q[1];
sx q[1];
rz(-0.97472753) q[1];
x q[2];
rz(-0.69334778) q[3];
sx q[3];
rz(-2.7487488) q[3];
sx q[3];
rz(-1.0885914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.990443) q[2];
sx q[2];
rz(-0.98439011) q[2];
sx q[2];
rz(-1.7669862) q[2];
rz(0.19045842) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16354887) q[0];
sx q[0];
rz(-1.6885641) q[0];
sx q[0];
rz(1.8769886) q[0];
rz(3.0679852) q[1];
sx q[1];
rz(-2.059748) q[1];
sx q[1];
rz(2.5216865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0479931) q[0];
sx q[0];
rz(-1.7353829) q[0];
sx q[0];
rz(-2.3175146) q[0];
x q[1];
rz(1.0629547) q[2];
sx q[2];
rz(-1.0217102) q[2];
sx q[2];
rz(-0.49875235) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.010308525) q[1];
sx q[1];
rz(-2.8493053) q[1];
sx q[1];
rz(0.75914219) q[1];
x q[2];
rz(-2.6750071) q[3];
sx q[3];
rz(-1.7853201) q[3];
sx q[3];
rz(1.9746575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2747043) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(0.14599027) q[2];
rz(-0.22249666) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(2.4116624) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20919007) q[0];
sx q[0];
rz(-1.1840191) q[0];
sx q[0];
rz(-2.9970072) q[0];
rz(1.9119541) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(0.59496224) q[2];
sx q[2];
rz(-1.7528201) q[2];
sx q[2];
rz(2.5217944) q[2];
rz(1.4208117) q[3];
sx q[3];
rz(-2.9290028) q[3];
sx q[3];
rz(1.3976189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
