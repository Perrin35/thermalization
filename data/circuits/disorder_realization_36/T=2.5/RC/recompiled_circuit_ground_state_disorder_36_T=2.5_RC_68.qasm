OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26389709) q[0];
sx q[0];
rz(-1.0385624) q[0];
sx q[0];
rz(-0.39499083) q[0];
rz(-2.360193) q[1];
sx q[1];
rz(-1.3716939) q[1];
sx q[1];
rz(-2.3514907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27777729) q[0];
sx q[0];
rz(-2.565747) q[0];
sx q[0];
rz(3.083549) q[0];
rz(-pi) q[1];
rz(-0.4460046) q[2];
sx q[2];
rz(-2.3141907) q[2];
sx q[2];
rz(2.3349886) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1829738) q[1];
sx q[1];
rz(-1.0369025) q[1];
sx q[1];
rz(0.16927232) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50689583) q[3];
sx q[3];
rz(-1.2032751) q[3];
sx q[3];
rz(1.3669922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6140952) q[2];
sx q[2];
rz(-3.0588394) q[2];
sx q[2];
rz(-0.56438524) q[2];
rz(2.2218521) q[3];
sx q[3];
rz(-0.65110937) q[3];
sx q[3];
rz(-1.1434327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1564002) q[0];
sx q[0];
rz(-0.99034482) q[0];
sx q[0];
rz(1.1457957) q[0];
rz(-1.0076373) q[1];
sx q[1];
rz(-0.8877019) q[1];
sx q[1];
rz(0.064858286) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.295395) q[0];
sx q[0];
rz(-1.0480651) q[0];
sx q[0];
rz(0.46434648) q[0];
rz(-pi) q[1];
rz(2.1880661) q[2];
sx q[2];
rz(-0.89260403) q[2];
sx q[2];
rz(1.8770521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8915598) q[1];
sx q[1];
rz(-1.3866826) q[1];
sx q[1];
rz(1.5450983) q[1];
rz(0.50527827) q[3];
sx q[3];
rz(-2.0483096) q[3];
sx q[3];
rz(-1.2569515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2603944) q[2];
sx q[2];
rz(-3.03646) q[2];
sx q[2];
rz(-0.51327389) q[2];
rz(-1.5567635) q[3];
sx q[3];
rz(-1.1791752) q[3];
sx q[3];
rz(-1.56196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2960812) q[0];
sx q[0];
rz(-0.054231461) q[0];
sx q[0];
rz(-0.84501141) q[0];
rz(0.71146479) q[1];
sx q[1];
rz(-0.80159801) q[1];
sx q[1];
rz(2.5149288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7539106) q[0];
sx q[0];
rz(-1.57608) q[0];
sx q[0];
rz(1.374701) q[0];
rz(-pi) q[1];
rz(1.1859796) q[2];
sx q[2];
rz(-1.3582412) q[2];
sx q[2];
rz(2.7351968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9756445) q[1];
sx q[1];
rz(-0.86417809) q[1];
sx q[1];
rz(2.9099275) q[1];
x q[2];
rz(-3.1342642) q[3];
sx q[3];
rz(-2.0617319) q[3];
sx q[3];
rz(2.1027264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6381548) q[2];
sx q[2];
rz(-1.2388836) q[2];
sx q[2];
rz(-2.9729291) q[2];
rz(-2.9851798) q[3];
sx q[3];
rz(-2.5037239) q[3];
sx q[3];
rz(0.28353459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6969358) q[0];
sx q[0];
rz(-2.038027) q[0];
sx q[0];
rz(2.7705833) q[0];
rz(-1.3320529) q[1];
sx q[1];
rz(-1.5947554) q[1];
sx q[1];
rz(2.7807049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7545032) q[0];
sx q[0];
rz(-1.4851804) q[0];
sx q[0];
rz(1.6667009) q[0];
rz(-pi) q[1];
rz(-1.5547916) q[2];
sx q[2];
rz(-1.4487378) q[2];
sx q[2];
rz(-2.6467488) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8239221) q[1];
sx q[1];
rz(-2.5054563) q[1];
sx q[1];
rz(-2.0348861) q[1];
rz(-pi) q[2];
rz(-2.4982483) q[3];
sx q[3];
rz(-0.78668252) q[3];
sx q[3];
rz(0.96970448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0648301) q[2];
sx q[2];
rz(-2.7702489) q[2];
sx q[2];
rz(1.4999464) q[2];
rz(-2.0867945) q[3];
sx q[3];
rz(-1.0902371) q[3];
sx q[3];
rz(1.1928026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37831369) q[0];
sx q[0];
rz(-1.167647) q[0];
sx q[0];
rz(-1.8462697) q[0];
rz(-0.12652215) q[1];
sx q[1];
rz(-0.69252068) q[1];
sx q[1];
rz(-1.893938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0678789) q[0];
sx q[0];
rz(-2.824221) q[0];
sx q[0];
rz(-2.6614891) q[0];
x q[1];
rz(-0.7681925) q[2];
sx q[2];
rz(-1.2821226) q[2];
sx q[2];
rz(1.4216258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8196053) q[1];
sx q[1];
rz(-1.7778559) q[1];
sx q[1];
rz(-2.7586069) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8528102) q[3];
sx q[3];
rz(-1.7502893) q[3];
sx q[3];
rz(3.0075551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7724472) q[2];
sx q[2];
rz(-0.36448604) q[2];
sx q[2];
rz(2.6920953) q[2];
rz(0.43826023) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(-1.0615758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19675572) q[0];
sx q[0];
rz(-2.0773092) q[0];
sx q[0];
rz(-2.7979895) q[0];
rz(-2.4876439) q[1];
sx q[1];
rz(-2.3908354) q[1];
sx q[1];
rz(0.59069815) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66309975) q[0];
sx q[0];
rz(-1.0859153) q[0];
sx q[0];
rz(-0.25867489) q[0];
rz(-2.3542064) q[2];
sx q[2];
rz(-2.2287268) q[2];
sx q[2];
rz(1.3546561) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6178083) q[1];
sx q[1];
rz(-1.7016546) q[1];
sx q[1];
rz(3.1310358) q[1];
rz(0.69575255) q[3];
sx q[3];
rz(-1.0067471) q[3];
sx q[3];
rz(-2.6524909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58064738) q[2];
sx q[2];
rz(-0.86981589) q[2];
sx q[2];
rz(-0.84442863) q[2];
rz(-1.0405509) q[3];
sx q[3];
rz(-1.2795762) q[3];
sx q[3];
rz(-1.9677264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3867253) q[0];
sx q[0];
rz(-2.6193021) q[0];
sx q[0];
rz(2.0544384) q[0];
rz(-2.6548751) q[1];
sx q[1];
rz(-1.646128) q[1];
sx q[1];
rz(-1.4901935) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9947858) q[0];
sx q[0];
rz(-1.675696) q[0];
sx q[0];
rz(-0.86284967) q[0];
x q[1];
rz(-1.0045687) q[2];
sx q[2];
rz(-0.90551584) q[2];
sx q[2];
rz(0.69349223) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5512276) q[1];
sx q[1];
rz(-1.1572954) q[1];
sx q[1];
rz(-1.1567924) q[1];
rz(-pi) q[2];
rz(-1.2734866) q[3];
sx q[3];
rz(-1.372223) q[3];
sx q[3];
rz(-1.8383593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7955486) q[2];
sx q[2];
rz(-0.86980021) q[2];
sx q[2];
rz(-2.8022433) q[2];
rz(-0.36335534) q[3];
sx q[3];
rz(-2.8715869) q[3];
sx q[3];
rz(1.5800765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5245847) q[0];
sx q[0];
rz(-2.0211077) q[0];
sx q[0];
rz(-2.8560915) q[0];
rz(-0.78954804) q[1];
sx q[1];
rz(-1.6098846) q[1];
sx q[1];
rz(-1.550536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11131903) q[0];
sx q[0];
rz(-0.34967129) q[0];
sx q[0];
rz(0.5446148) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0976426) q[2];
sx q[2];
rz(-1.6990597) q[2];
sx q[2];
rz(-1.8086901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4553677) q[1];
sx q[1];
rz(-1.1535201) q[1];
sx q[1];
rz(-2.5646006) q[1];
rz(-pi) q[2];
rz(0.66724369) q[3];
sx q[3];
rz(-1.5074456) q[3];
sx q[3];
rz(-2.1433307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8153814) q[2];
sx q[2];
rz(-2.9823494) q[2];
sx q[2];
rz(-0.86686575) q[2];
rz(0.73831144) q[3];
sx q[3];
rz(-1.5840931) q[3];
sx q[3];
rz(0.90832925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065780491) q[0];
sx q[0];
rz(-0.74956885) q[0];
sx q[0];
rz(-0.69920364) q[0];
rz(2.5993644) q[1];
sx q[1];
rz(-1.3424073) q[1];
sx q[1];
rz(2.8453765) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4980221) q[0];
sx q[0];
rz(-1.5012465) q[0];
sx q[0];
rz(-0.83840094) q[0];
x q[1];
rz(2.3260957) q[2];
sx q[2];
rz(-1.0692714) q[2];
sx q[2];
rz(3.0072711) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6449127) q[1];
sx q[1];
rz(-1.0154004) q[1];
sx q[1];
rz(2.4128561) q[1];
x q[2];
rz(1.1929729) q[3];
sx q[3];
rz(-1.4903896) q[3];
sx q[3];
rz(-0.37624826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88973033) q[2];
sx q[2];
rz(-2.6984062) q[2];
sx q[2];
rz(2.6207793) q[2];
rz(-0.03171799) q[3];
sx q[3];
rz(-2.7364276) q[3];
sx q[3];
rz(-3.090455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6601335) q[0];
sx q[0];
rz(-2.0429459) q[0];
sx q[0];
rz(-2.3302186) q[0];
rz(-0.42586455) q[1];
sx q[1];
rz(-2.2399797) q[1];
sx q[1];
rz(-1.834555) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15430121) q[0];
sx q[0];
rz(-1.5736921) q[0];
sx q[0];
rz(1.5737246) q[0];
rz(-pi) q[1];
rz(1.6128236) q[2];
sx q[2];
rz(-2.6767455) q[2];
sx q[2];
rz(2.5744409) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.934224) q[1];
sx q[1];
rz(-1.1282053) q[1];
sx q[1];
rz(1.9277444) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7440935) q[3];
sx q[3];
rz(-1.8911282) q[3];
sx q[3];
rz(-0.082435247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9797152) q[2];
sx q[2];
rz(-2.3780509) q[2];
sx q[2];
rz(2.7757576) q[2];
rz(2.2063935) q[3];
sx q[3];
rz(-1.5949944) q[3];
sx q[3];
rz(-0.38177761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314421) q[0];
sx q[0];
rz(-1.3642385) q[0];
sx q[0];
rz(1.4862899) q[0];
rz(-1.4797795) q[1];
sx q[1];
rz(-1.9098837) q[1];
sx q[1];
rz(2.2178537) q[1];
rz(-2.1187028) q[2];
sx q[2];
rz(-1.6216424) q[2];
sx q[2];
rz(-2.4301823) q[2];
rz(1.4167841) q[3];
sx q[3];
rz(-1.8324251) q[3];
sx q[3];
rz(0.13941924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
