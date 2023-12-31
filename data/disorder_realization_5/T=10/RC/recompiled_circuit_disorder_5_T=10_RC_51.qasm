OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(4.7441626) q[0];
sx q[0];
rz(6.8466853) q[0];
sx q[0];
rz(9.8817649) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5029992) q[0];
sx q[0];
rz(-2.8832957) q[0];
sx q[0];
rz(0.55708416) q[0];
x q[1];
rz(1.4957501) q[2];
sx q[2];
rz(-1.6610166) q[2];
sx q[2];
rz(0.98352945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(-3.0800746) q[1];
x q[2];
rz(0.92313487) q[3];
sx q[3];
rz(-1.2520408) q[3];
sx q[3];
rz(-1.8766599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-0.045036137) q[0];
rz(2.8024407) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4883603) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(-1.0351719) q[0];
rz(-pi) q[1];
rz(0.10404189) q[2];
sx q[2];
rz(-1.0696971) q[2];
sx q[2];
rz(0.83545557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3992577) q[1];
sx q[1];
rz(-2.5207673) q[1];
sx q[1];
rz(1.6595112) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0117399) q[3];
sx q[3];
rz(-2.2359071) q[3];
sx q[3];
rz(-3.0761884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(-1.7155898) q[2];
rz(0.60570335) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.3396938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08683603) q[0];
sx q[0];
rz(-1.1491547) q[0];
sx q[0];
rz(-0.7936759) q[0];
x q[1];
rz(0.3091829) q[2];
sx q[2];
rz(-0.90484607) q[2];
sx q[2];
rz(3.0004629) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(-1.3819441) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6684746) q[3];
sx q[3];
rz(-2.3948673) q[3];
sx q[3];
rz(-2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(2.9273422) q[0];
rz(-2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(-2.4750211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4857805) q[0];
sx q[0];
rz(-2.862145) q[0];
sx q[0];
rz(-0.73894545) q[0];
rz(-pi) q[1];
rz(-2.8104086) q[2];
sx q[2];
rz(-0.74197717) q[2];
sx q[2];
rz(-0.77510288) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9208593) q[1];
sx q[1];
rz(-1.024106) q[1];
sx q[1];
rz(2.2842555) q[1];
rz(-pi) q[2];
rz(0.14931071) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(-3.1130476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(0.66037035) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(2.3045325) q[0];
rz(2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.1594835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75393049) q[0];
sx q[0];
rz(-1.2571063) q[0];
sx q[0];
rz(2.4051106) q[0];
rz(-0.97255623) q[2];
sx q[2];
rz(-1.4417357) q[2];
sx q[2];
rz(-1.8051403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24134484) q[1];
sx q[1];
rz(-0.7175788) q[1];
sx q[1];
rz(-2.248583) q[1];
rz(2.106835) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-3.0767373) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(2.0264453) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(-3.040722) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3374702) q[0];
sx q[0];
rz(-0.94479783) q[0];
sx q[0];
rz(-1.3834329) q[0];
rz(-pi) q[1];
rz(-1.8632554) q[2];
sx q[2];
rz(-2.0403701) q[2];
sx q[2];
rz(-1.3172319) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3979891) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(2.9823751) q[1];
rz(-0.40525873) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(-2.1586771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(-1.8484176) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9822385) q[0];
sx q[0];
rz(-1.9040477) q[0];
sx q[0];
rz(-0.23157816) q[0];
rz(-pi) q[1];
rz(-2.1805448) q[2];
sx q[2];
rz(-0.38182048) q[2];
sx q[2];
rz(1.0605937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75377611) q[1];
sx q[1];
rz(-0.13145914) q[1];
sx q[1];
rz(-0.17863518) q[1];
rz(-1.5457821) q[3];
sx q[3];
rz(-0.56846148) q[3];
sx q[3];
rz(-0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-3.0623073) q[0];
rz(-1.1212564) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.1987196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2326956) q[0];
sx q[0];
rz(-2.7017653) q[0];
sx q[0];
rz(-2.9798685) q[0];
rz(-pi) q[1];
rz(1.347581) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(-1.5471293) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.63003507) q[1];
sx q[1];
rz(-2.6361598) q[1];
sx q[1];
rz(3.0384484) q[1];
rz(-pi) q[2];
rz(-2.6783887) q[3];
sx q[3];
rz(-1.8052881) q[3];
sx q[3];
rz(0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85598677) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-0.58832204) q[0];
rz(-3.1148124) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(2.2081597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8342455) q[0];
sx q[0];
rz(-1.5114307) q[0];
sx q[0];
rz(-0.10660118) q[0];
rz(0.33414267) q[2];
sx q[2];
rz(-1.9311937) q[2];
sx q[2];
rz(0.44597086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48298207) q[1];
sx q[1];
rz(-0.88739363) q[1];
sx q[1];
rz(-1.618209) q[1];
x q[2];
rz(0.65255717) q[3];
sx q[3];
rz(-1.2823294) q[3];
sx q[3];
rz(-2.1289338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.062722) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(0.36994568) q[2];
rz(2.2423819) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(-1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-2.9577589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18098022) q[0];
sx q[0];
rz(-1.6913337) q[0];
sx q[0];
rz(-1.6880077) q[0];
rz(-pi) q[1];
rz(-2.3046938) q[2];
sx q[2];
rz(-2.1914334) q[2];
sx q[2];
rz(-1.2037954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.373133) q[1];
sx q[1];
rz(-0.8994973) q[1];
sx q[1];
rz(-0.21233227) q[1];
x q[2];
rz(0.023601836) q[3];
sx q[3];
rz(-0.58348237) q[3];
sx q[3];
rz(-0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(0.72485483) q[2];
sx q[2];
rz(-0.58314322) q[2];
sx q[2];
rz(-2.4287139) q[2];
rz(-1.3115053) q[3];
sx q[3];
rz(-0.6060096) q[3];
sx q[3];
rz(0.93939645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
