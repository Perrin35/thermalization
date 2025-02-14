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
rz(0.67233664) q[0];
sx q[0];
rz(5.0007102) q[0];
sx q[0];
rz(11.034427) q[0];
rz(-1.9077644) q[1];
sx q[1];
rz(-1.2169853) q[1];
sx q[1];
rz(1.4438862) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1506526) q[0];
sx q[0];
rz(-0.67545891) q[0];
sx q[0];
rz(1.7016241) q[0];
rz(0.89272372) q[2];
sx q[2];
rz(-0.44995445) q[2];
sx q[2];
rz(-1.3154678) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68221625) q[1];
sx q[1];
rz(-0.77118783) q[1];
sx q[1];
rz(-0.60318635) q[1];
x q[2];
rz(1.4920294) q[3];
sx q[3];
rz(-0.34011671) q[3];
sx q[3];
rz(2.801737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7794789) q[2];
sx q[2];
rz(-1.2059836) q[2];
sx q[2];
rz(-2.3625679) q[2];
rz(-1.3075167) q[3];
sx q[3];
rz(-0.3236168) q[3];
sx q[3];
rz(-1.754508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5883412) q[0];
sx q[0];
rz(-1.3655751) q[0];
sx q[0];
rz(-2.8643082) q[0];
rz(-2.7514027) q[1];
sx q[1];
rz(-1.1017825) q[1];
sx q[1];
rz(2.0139093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41306333) q[0];
sx q[0];
rz(-1.5507694) q[0];
sx q[0];
rz(-0.042547556) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28553005) q[2];
sx q[2];
rz(-1.4737616) q[2];
sx q[2];
rz(2.3842528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0223282) q[1];
sx q[1];
rz(-2.8624967) q[1];
sx q[1];
rz(-1.7995681) q[1];
rz(-0.47445837) q[3];
sx q[3];
rz(-2.4224307) q[3];
sx q[3];
rz(-0.31960318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7057425) q[2];
sx q[2];
rz(-2.5472992) q[2];
sx q[2];
rz(1.1951813) q[2];
rz(2.4863906) q[3];
sx q[3];
rz(-1.4846669) q[3];
sx q[3];
rz(-0.5932194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6708267) q[0];
sx q[0];
rz(-0.68597454) q[0];
sx q[0];
rz(2.2210448) q[0];
rz(1.2819194) q[1];
sx q[1];
rz(-1.1290519) q[1];
sx q[1];
rz(-1.4528073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6920928) q[0];
sx q[0];
rz(-1.458174) q[0];
sx q[0];
rz(-1.4829259) q[0];
rz(0.31450504) q[2];
sx q[2];
rz(-2.1172519) q[2];
sx q[2];
rz(-0.98775452) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9640426) q[1];
sx q[1];
rz(-1.519555) q[1];
sx q[1];
rz(-2.8058626) q[1];
rz(-0.7887855) q[3];
sx q[3];
rz(-0.093063912) q[3];
sx q[3];
rz(-2.399202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61536962) q[2];
sx q[2];
rz(-2.0873895) q[2];
sx q[2];
rz(-0.39895454) q[2];
rz(-2.6279972) q[3];
sx q[3];
rz(-2.8245638) q[3];
sx q[3];
rz(-1.1494466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8154163) q[0];
sx q[0];
rz(-2.5626817) q[0];
sx q[0];
rz(1.2166566) q[0];
rz(-2.6293964) q[1];
sx q[1];
rz(-2.0307816) q[1];
sx q[1];
rz(-1.1042575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51795635) q[0];
sx q[0];
rz(-0.41455634) q[0];
sx q[0];
rz(-0.74409938) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7521367) q[2];
sx q[2];
rz(-1.7415819) q[2];
sx q[2];
rz(2.6362399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0287289) q[1];
sx q[1];
rz(-2.9359438) q[1];
sx q[1];
rz(-1.2713232) q[1];
rz(-pi) q[2];
rz(2.9151408) q[3];
sx q[3];
rz(-1.8237916) q[3];
sx q[3];
rz(-1.7915726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5451374) q[2];
sx q[2];
rz(-1.9497654) q[2];
sx q[2];
rz(-2.5917501) q[2];
rz(0.95692974) q[3];
sx q[3];
rz(-2.3423829) q[3];
sx q[3];
rz(1.6526615) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6560219) q[0];
sx q[0];
rz(-0.85000426) q[0];
sx q[0];
rz(-0.75918424) q[0];
rz(-2.1929269) q[1];
sx q[1];
rz(-1.1439088) q[1];
sx q[1];
rz(-2.8099828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.054197) q[0];
sx q[0];
rz(-2.3080462) q[0];
sx q[0];
rz(1.1655318) q[0];
rz(-0.070721432) q[2];
sx q[2];
rz(-0.47571936) q[2];
sx q[2];
rz(-2.4782319) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4840821) q[1];
sx q[1];
rz(-0.23536029) q[1];
sx q[1];
rz(2.2886307) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7071574) q[3];
sx q[3];
rz(-1.9788187) q[3];
sx q[3];
rz(1.7849762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18629508) q[2];
sx q[2];
rz(-2.1897327) q[2];
sx q[2];
rz(2.1985506) q[2];
rz(-2.9694563) q[3];
sx q[3];
rz(-0.94477263) q[3];
sx q[3];
rz(-0.22323639) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73582369) q[0];
sx q[0];
rz(-0.74546927) q[0];
sx q[0];
rz(-1.5663585) q[0];
rz(2.9834421) q[1];
sx q[1];
rz(-0.76845303) q[1];
sx q[1];
rz(0.92811531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0029308) q[0];
sx q[0];
rz(-1.2671906) q[0];
sx q[0];
rz(-1.7055442) q[0];
x q[1];
rz(2.1238223) q[2];
sx q[2];
rz(-1.206351) q[2];
sx q[2];
rz(1.6318969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21213046) q[1];
sx q[1];
rz(-1.3155335) q[1];
sx q[1];
rz(-0.90025951) q[1];
x q[2];
rz(-0.35291725) q[3];
sx q[3];
rz(-1.4541876) q[3];
sx q[3];
rz(1.1427134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9485665) q[2];
sx q[2];
rz(-2.8552449) q[2];
sx q[2];
rz(2.7890653) q[2];
rz(2.7726717) q[3];
sx q[3];
rz(-0.58266321) q[3];
sx q[3];
rz(0.85744706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3567268) q[0];
sx q[0];
rz(-3.1251188) q[0];
sx q[0];
rz(0.67100588) q[0];
rz(0.10637936) q[1];
sx q[1];
rz(-2.2439067) q[1];
sx q[1];
rz(0.62972128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8917568) q[0];
sx q[0];
rz(-1.5895004) q[0];
sx q[0];
rz(2.9500701) q[0];
x q[1];
rz(-0.81083576) q[2];
sx q[2];
rz(-0.5631643) q[2];
sx q[2];
rz(0.25897464) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32618648) q[1];
sx q[1];
rz(-1.2478653) q[1];
sx q[1];
rz(-1.4923444) q[1];
rz(-pi) q[2];
rz(2.4829818) q[3];
sx q[3];
rz(-1.8192184) q[3];
sx q[3];
rz(-0.90693865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.070039198) q[2];
sx q[2];
rz(-0.89954251) q[2];
sx q[2];
rz(2.1001935) q[2];
rz(-1.533016) q[3];
sx q[3];
rz(-0.93419111) q[3];
sx q[3];
rz(-1.9937203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94569412) q[0];
sx q[0];
rz(-2.4905289) q[0];
sx q[0];
rz(-2.661929) q[0];
rz(0.02462968) q[1];
sx q[1];
rz(-1.004091) q[1];
sx q[1];
rz(3.0020795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280645) q[0];
sx q[0];
rz(-1.6798927) q[0];
sx q[0];
rz(0.086221545) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37711511) q[2];
sx q[2];
rz(-2.0465436) q[2];
sx q[2];
rz(-1.7390133) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3132726) q[1];
sx q[1];
rz(-1.4045776) q[1];
sx q[1];
rz(0.66201026) q[1];
rz(-pi) q[2];
rz(0.068753763) q[3];
sx q[3];
rz(-0.67151755) q[3];
sx q[3];
rz(0.17291394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2554243) q[2];
sx q[2];
rz(-1.5972127) q[2];
sx q[2];
rz(1.9891948) q[2];
rz(1.7216916) q[3];
sx q[3];
rz(-0.73437381) q[3];
sx q[3];
rz(-1.3091807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060147978) q[0];
sx q[0];
rz(-1.4652493) q[0];
sx q[0];
rz(1.4981221) q[0];
rz(2.3807047) q[1];
sx q[1];
rz(-0.94327578) q[1];
sx q[1];
rz(1.6852185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0022572) q[0];
sx q[0];
rz(-0.59334785) q[0];
sx q[0];
rz(-0.32440455) q[0];
rz(-pi) q[1];
rz(-2.7414906) q[2];
sx q[2];
rz(-1.3753554) q[2];
sx q[2];
rz(0.2646499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4059116) q[1];
sx q[1];
rz(-1.0764437) q[1];
sx q[1];
rz(0.45208652) q[1];
x q[2];
rz(-0.96181576) q[3];
sx q[3];
rz(-2.4815593) q[3];
sx q[3];
rz(-0.48840085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8196408) q[2];
sx q[2];
rz(-2.7072622) q[2];
sx q[2];
rz(1.0158018) q[2];
rz(2.2470233) q[3];
sx q[3];
rz(-1.8958478) q[3];
sx q[3];
rz(2.3999124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9858518) q[0];
sx q[0];
rz(-2.4116801) q[0];
sx q[0];
rz(2.1515382) q[0];
rz(-2.2186225) q[1];
sx q[1];
rz(-1.3786517) q[1];
sx q[1];
rz(2.1754481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7434563) q[0];
sx q[0];
rz(-1.4706637) q[0];
sx q[0];
rz(0.070835367) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82426333) q[2];
sx q[2];
rz(-2.5474862) q[2];
sx q[2];
rz(2.4115393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6412539) q[1];
sx q[1];
rz(-2.6473844) q[1];
sx q[1];
rz(0.92730913) q[1];
rz(1.0564141) q[3];
sx q[3];
rz(-2.9136411) q[3];
sx q[3];
rz(-2.0755656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82603377) q[2];
sx q[2];
rz(-0.78353271) q[2];
sx q[2];
rz(-2.3890736) q[2];
rz(-3.0044532) q[3];
sx q[3];
rz(-0.64118853) q[3];
sx q[3];
rz(0.34998676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.8990477) q[0];
sx q[0];
rz(-2.516863) q[0];
sx q[0];
rz(2.94009) q[0];
rz(-1.7789727) q[1];
sx q[1];
rz(-2.2687804) q[1];
sx q[1];
rz(-0.16558095) q[1];
rz(3.0918783) q[2];
sx q[2];
rz(-1.4820953) q[2];
sx q[2];
rz(-2.6470071) q[2];
rz(2.2929819) q[3];
sx q[3];
rz(-0.52042421) q[3];
sx q[3];
rz(2.0757794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
