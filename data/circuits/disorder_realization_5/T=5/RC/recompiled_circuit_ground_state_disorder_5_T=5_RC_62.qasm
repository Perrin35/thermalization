OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(-2.2192686) q[0];
rz(-0.84696472) q[1];
sx q[1];
rz(4.6159336) q[1];
sx q[1];
rz(9.6428975) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51061741) q[0];
sx q[0];
rz(-1.4790863) q[0];
sx q[0];
rz(-2.6789078) q[0];
rz(-pi) q[1];
rz(-1.6118902) q[2];
sx q[2];
rz(-1.9485222) q[2];
sx q[2];
rz(2.1884577) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.328124) q[1];
sx q[1];
rz(-1.3513997) q[1];
sx q[1];
rz(0.043299999) q[1];
rz(-pi) q[2];
rz(-1.4014729) q[3];
sx q[3];
rz(-2.7062731) q[3];
sx q[3];
rz(3.0141413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0685136) q[2];
sx q[2];
rz(-1.2129236) q[2];
sx q[2];
rz(2.6050513) q[2];
rz(2.1327298) q[3];
sx q[3];
rz(-2.3705685) q[3];
sx q[3];
rz(-2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63758481) q[0];
sx q[0];
rz(-1.3115839) q[0];
sx q[0];
rz(-3.1245533) q[0];
rz(1.0631961) q[1];
sx q[1];
rz(-0.76779643) q[1];
sx q[1];
rz(-0.95471901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9794036) q[0];
sx q[0];
rz(-1.6260176) q[0];
sx q[0];
rz(-0.22947854) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9770697) q[2];
sx q[2];
rz(-1.2753107) q[2];
sx q[2];
rz(-2.4128259) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8361289) q[1];
sx q[1];
rz(-1.8082128) q[1];
sx q[1];
rz(1.8029455) q[1];
rz(-pi) q[2];
rz(1.6240261) q[3];
sx q[3];
rz(-1.551515) q[3];
sx q[3];
rz(-2.791491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9118328) q[2];
sx q[2];
rz(-0.91270295) q[2];
sx q[2];
rz(2.0274053) q[2];
rz(2.0071425) q[3];
sx q[3];
rz(-0.19616923) q[3];
sx q[3];
rz(0.1964868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5512307) q[0];
sx q[0];
rz(-0.95142618) q[0];
sx q[0];
rz(0.1828585) q[0];
rz(-3.0602449) q[1];
sx q[1];
rz(-2.3902939) q[1];
sx q[1];
rz(0.61838165) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8066766) q[0];
sx q[0];
rz(-0.78640079) q[0];
sx q[0];
rz(-1.4674835) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3803394) q[2];
sx q[2];
rz(-2.3229685) q[2];
sx q[2];
rz(-3.034575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.394819) q[1];
sx q[1];
rz(-1.9095943) q[1];
sx q[1];
rz(2.9940786) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86966536) q[3];
sx q[3];
rz(-2.0732862) q[3];
sx q[3];
rz(-1.089407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0028093) q[2];
sx q[2];
rz(-1.8178136) q[2];
sx q[2];
rz(2.7749824) q[2];
rz(-1.2159411) q[3];
sx q[3];
rz(-1.1953019) q[3];
sx q[3];
rz(0.6634357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65130305) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(-2.41462) q[0];
rz(-2.2333249) q[1];
sx q[1];
rz(-0.5779225) q[1];
sx q[1];
rz(-1.04331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3573049) q[0];
sx q[0];
rz(-2.4803964) q[0];
sx q[0];
rz(0.0088328514) q[0];
rz(-pi) q[1];
rz(-0.06177549) q[2];
sx q[2];
rz(-1.1959658) q[2];
sx q[2];
rz(1.0967364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.21884313) q[1];
sx q[1];
rz(-2.6237539) q[1];
sx q[1];
rz(-1.0475558) q[1];
x q[2];
rz(0.38049145) q[3];
sx q[3];
rz(-2.0999319) q[3];
sx q[3];
rz(0.38705119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72994453) q[2];
sx q[2];
rz(-2.8204155) q[2];
sx q[2];
rz(-1.3058861) q[2];
rz(-0.11792396) q[3];
sx q[3];
rz(-2.6730461) q[3];
sx q[3];
rz(2.0809295) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0333198) q[0];
sx q[0];
rz(-2.1554027) q[0];
sx q[0];
rz(-1.5340075) q[0];
rz(-2.8458505) q[1];
sx q[1];
rz(-1.60138) q[1];
sx q[1];
rz(-1.6827513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47285493) q[0];
sx q[0];
rz(-0.74873996) q[0];
sx q[0];
rz(1.4677901) q[0];
rz(-pi) q[1];
rz(2.8334728) q[2];
sx q[2];
rz(-1.9240555) q[2];
sx q[2];
rz(2.0120914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.69469317) q[1];
sx q[1];
rz(-2.2647479) q[1];
sx q[1];
rz(2.8847221) q[1];
rz(-pi) q[2];
rz(-1.9146054) q[3];
sx q[3];
rz(-0.60887486) q[3];
sx q[3];
rz(2.6725149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7532588) q[2];
sx q[2];
rz(-2.7497141) q[2];
sx q[2];
rz(2.4957116) q[2];
rz(-1.9258457) q[3];
sx q[3];
rz(-1.682621) q[3];
sx q[3];
rz(1.7997883) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011768613) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(2.5339793) q[0];
rz(1.1786849) q[1];
sx q[1];
rz(-1.2858398) q[1];
sx q[1];
rz(-2.7778621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54485142) q[0];
sx q[0];
rz(-2.5646696) q[0];
sx q[0];
rz(-2.0175009) q[0];
rz(-0.046437101) q[2];
sx q[2];
rz(-1.9032459) q[2];
sx q[2];
rz(-0.22975555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44139578) q[1];
sx q[1];
rz(-1.6906941) q[1];
sx q[1];
rz(-2.3223206) q[1];
rz(-pi) q[2];
x q[2];
rz(2.167114) q[3];
sx q[3];
rz(-1.5502366) q[3];
sx q[3];
rz(-2.1316776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5522449) q[2];
sx q[2];
rz(-3.037368) q[2];
sx q[2];
rz(-1.1927401) q[2];
rz(-0.73181152) q[3];
sx q[3];
rz(-1.4251499) q[3];
sx q[3];
rz(1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.28855395) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(1.5227675) q[0];
rz(1.6743926) q[1];
sx q[1];
rz(-1.3102945) q[1];
sx q[1];
rz(2.4640962) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38792363) q[0];
sx q[0];
rz(-0.69783995) q[0];
sx q[0];
rz(1.4927255) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4676553) q[2];
sx q[2];
rz(-0.29602414) q[2];
sx q[2];
rz(0.37002555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4851393) q[1];
sx q[1];
rz(-2.5514126) q[1];
sx q[1];
rz(0.82832576) q[1];
rz(1.4114289) q[3];
sx q[3];
rz(-0.84468088) q[3];
sx q[3];
rz(2.1350525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39685321) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(1.8358561) q[2];
rz(-0.33852494) q[3];
sx q[3];
rz(-2.0207696) q[3];
sx q[3];
rz(-0.6849851) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5386388) q[0];
sx q[0];
rz(-2.9260577) q[0];
sx q[0];
rz(2.9576874) q[0];
rz(-1.454608) q[1];
sx q[1];
rz(-2.589476) q[1];
sx q[1];
rz(-2.6457381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5372779) q[0];
sx q[0];
rz(-1.0838944) q[0];
sx q[0];
rz(-2.8394798) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3681201) q[2];
sx q[2];
rz(-2.2272553) q[2];
sx q[2];
rz(0.17420775) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4019805) q[1];
sx q[1];
rz(-2.6079) q[1];
sx q[1];
rz(2.3237565) q[1];
rz(-pi) q[2];
rz(0.11885507) q[3];
sx q[3];
rz(-2.1200075) q[3];
sx q[3];
rz(2.1191747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0869861) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(1.1620713) q[2];
rz(1.6880796) q[3];
sx q[3];
rz(-0.96265692) q[3];
sx q[3];
rz(1.8614205) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7206955) q[0];
sx q[0];
rz(-1.9891885) q[0];
sx q[0];
rz(-2.5757117) q[0];
rz(0.3903009) q[1];
sx q[1];
rz(-1.5696328) q[1];
sx q[1];
rz(-1.8338667) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6193251) q[0];
sx q[0];
rz(-0.59460708) q[0];
sx q[0];
rz(2.1451166) q[0];
x q[1];
rz(0.21251596) q[2];
sx q[2];
rz(-1.9409928) q[2];
sx q[2];
rz(-3.1345362) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3257287) q[1];
sx q[1];
rz(-2.1727409) q[1];
sx q[1];
rz(-1.25156) q[1];
rz(-pi) q[2];
rz(-0.89687895) q[3];
sx q[3];
rz(-0.19052902) q[3];
sx q[3];
rz(-2.6485557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86964837) q[2];
sx q[2];
rz(-0.41676909) q[2];
sx q[2];
rz(-1.0375674) q[2];
rz(1.8218254) q[3];
sx q[3];
rz(-2.3128553) q[3];
sx q[3];
rz(-2.493609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8808402) q[0];
sx q[0];
rz(-0.97761959) q[0];
sx q[0];
rz(-2.4993437) q[0];
rz(-1.3308659) q[1];
sx q[1];
rz(-2.2763177) q[1];
sx q[1];
rz(-0.16255249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560682) q[0];
sx q[0];
rz(-2.3352288) q[0];
sx q[0];
rz(-2.1245405) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3847975) q[2];
sx q[2];
rz(-2.6399603) q[2];
sx q[2];
rz(2.121814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.049920883) q[1];
sx q[1];
rz(-2.0494048) q[1];
sx q[1];
rz(0.19772999) q[1];
rz(-pi) q[2];
rz(2.558296) q[3];
sx q[3];
rz(-1.5044893) q[3];
sx q[3];
rz(2.5401859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51005298) q[2];
sx q[2];
rz(-1.2756462) q[2];
sx q[2];
rz(-2.1007382) q[2];
rz(2.1736274) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(2.4805243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(2.9087635) q[1];
sx q[1];
rz(-2.1131344) q[1];
sx q[1];
rz(0.0075385787) q[1];
rz(0.54388028) q[2];
sx q[2];
rz(-0.82533045) q[2];
sx q[2];
rz(-2.5173204) q[2];
rz(1.1875752) q[3];
sx q[3];
rz(-2.2338727) q[3];
sx q[3];
rz(3.0919872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
