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
rz(-3.0109316) q[0];
sx q[0];
rz(-1.8008512) q[0];
sx q[0];
rz(-1.2641579) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(-1.0882508) q[1];
sx q[1];
rz(0.01297125) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54161763) q[0];
sx q[0];
rz(-2.0429039) q[0];
sx q[0];
rz(-1.7469445) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58836909) q[2];
sx q[2];
rz(-2.9466963) q[2];
sx q[2];
rz(2.4590059) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4526635) q[1];
sx q[1];
rz(-2.6580992) q[1];
sx q[1];
rz(-0.74093282) q[1];
x q[2];
rz(2.8230328) q[3];
sx q[3];
rz(-2.5346642) q[3];
sx q[3];
rz(-1.7817117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26244792) q[2];
sx q[2];
rz(-1.6877354) q[2];
sx q[2];
rz(-3.0461779) q[2];
rz(-2.6573507) q[3];
sx q[3];
rz(-0.80317322) q[3];
sx q[3];
rz(1.4580844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478304) q[0];
sx q[0];
rz(-1.4365124) q[0];
sx q[0];
rz(-2.4875212) q[0];
rz(1.7895128) q[1];
sx q[1];
rz(-2.4857931) q[1];
sx q[1];
rz(-0.10993122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9655361) q[0];
sx q[0];
rz(-1.5829464) q[0];
sx q[0];
rz(-1.5909082) q[0];
rz(2.201033) q[2];
sx q[2];
rz(-1.1917243) q[2];
sx q[2];
rz(-2.4565168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0866894) q[1];
sx q[1];
rz(-0.48799054) q[1];
sx q[1];
rz(-2.7665374) q[1];
rz(-pi) q[2];
rz(1.3160922) q[3];
sx q[3];
rz(-1.1858018) q[3];
sx q[3];
rz(0.32190548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8249417) q[2];
sx q[2];
rz(-1.9330838) q[2];
sx q[2];
rz(-1.4671154) q[2];
rz(2.1064827) q[3];
sx q[3];
rz(-0.78101522) q[3];
sx q[3];
rz(-1.8234113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3159981) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(0.79063928) q[0];
rz(0.74360338) q[1];
sx q[1];
rz(-1.5433106) q[1];
sx q[1];
rz(0.57949439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7046708) q[0];
sx q[0];
rz(-1.0365067) q[0];
sx q[0];
rz(2.2493258) q[0];
rz(-1.2762464) q[2];
sx q[2];
rz(-1.8229457) q[2];
sx q[2];
rz(-1.1974826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0620561) q[1];
sx q[1];
rz(-1.3128451) q[1];
sx q[1];
rz(0.038611488) q[1];
x q[2];
rz(-2.4030645) q[3];
sx q[3];
rz(-0.66940847) q[3];
sx q[3];
rz(0.51046342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6774595) q[2];
sx q[2];
rz(-2.3049998) q[2];
sx q[2];
rz(0.38132384) q[2];
rz(2.2903806) q[3];
sx q[3];
rz(-0.92654735) q[3];
sx q[3];
rz(0.73494953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320936) q[0];
sx q[0];
rz(-2.5272326) q[0];
sx q[0];
rz(0.36439782) q[0];
rz(2.9413307) q[1];
sx q[1];
rz(-1.7297144) q[1];
sx q[1];
rz(-0.099954896) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1276949) q[0];
sx q[0];
rz(-1.5504) q[0];
sx q[0];
rz(1.4958032) q[0];
rz(-1.7968788) q[2];
sx q[2];
rz(-0.43856171) q[2];
sx q[2];
rz(-1.3815224) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.050592) q[1];
sx q[1];
rz(-3.0464131) q[1];
sx q[1];
rz(-0.23732136) q[1];
x q[2];
rz(-2.1853946) q[3];
sx q[3];
rz(-2.0609988) q[3];
sx q[3];
rz(-0.7008926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.072307) q[2];
sx q[2];
rz(-2.0733209) q[2];
sx q[2];
rz(-0.11492534) q[2];
rz(-2.0011486) q[3];
sx q[3];
rz(-1.2830257) q[3];
sx q[3];
rz(-2.1663402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.2045778) q[0];
sx q[0];
rz(-0.95252043) q[0];
sx q[0];
rz(-2.9691147) q[0];
rz(1.0109488) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(0.15484658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60733838) q[0];
sx q[0];
rz(-1.3528498) q[0];
sx q[0];
rz(1.2675257) q[0];
rz(0.24598083) q[2];
sx q[2];
rz(-2.2005251) q[2];
sx q[2];
rz(-0.069815947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46609572) q[1];
sx q[1];
rz(-1.4267529) q[1];
sx q[1];
rz(-0.41464582) q[1];
rz(-pi) q[2];
rz(1.6825292) q[3];
sx q[3];
rz(-1.7479241) q[3];
sx q[3];
rz(0.88645173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0871206) q[2];
sx q[2];
rz(-0.27721578) q[2];
sx q[2];
rz(2.7692914) q[2];
rz(-2.7436658) q[3];
sx q[3];
rz(-1.1436661) q[3];
sx q[3];
rz(0.00042375617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.953124) q[0];
sx q[0];
rz(-2.5690881) q[0];
sx q[0];
rz(1.5931801) q[0];
rz(-2.7717223) q[1];
sx q[1];
rz(-1.7431755) q[1];
sx q[1];
rz(-2.5028548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4700482) q[0];
sx q[0];
rz(-0.28451583) q[0];
sx q[0];
rz(-2.1485062) q[0];
rz(-2.0962786) q[2];
sx q[2];
rz(-1.3945701) q[2];
sx q[2];
rz(2.3228563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7228917) q[1];
sx q[1];
rz(-1.9564374) q[1];
sx q[1];
rz(-2.0608749) q[1];
rz(-pi) q[2];
x q[2];
rz(1.005976) q[3];
sx q[3];
rz(-2.2788439) q[3];
sx q[3];
rz(-2.4410409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0686191) q[2];
sx q[2];
rz(-2.0899453) q[2];
sx q[2];
rz(2.8886786) q[2];
rz(-2.2933293) q[3];
sx q[3];
rz(-1.4784644) q[3];
sx q[3];
rz(-3.0566791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0406168) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(2.9926391) q[0];
rz(1.8303998) q[1];
sx q[1];
rz(-1.7313749) q[1];
sx q[1];
rz(-0.054100903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9210756) q[0];
sx q[0];
rz(-2.069325) q[0];
sx q[0];
rz(-1.4709298) q[0];
rz(-pi) q[1];
rz(-3.0981423) q[2];
sx q[2];
rz(-2.3386526) q[2];
sx q[2];
rz(-0.96722764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3971567) q[1];
sx q[1];
rz(-0.87578339) q[1];
sx q[1];
rz(1.4276464) q[1];
rz(2.0460753) q[3];
sx q[3];
rz(-1.5500808) q[3];
sx q[3];
rz(1.4521862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3397303) q[2];
sx q[2];
rz(-1.9714377) q[2];
sx q[2];
rz(-0.96699634) q[2];
rz(0.70080924) q[3];
sx q[3];
rz(-2.1613224) q[3];
sx q[3];
rz(-0.35082671) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0130149) q[0];
sx q[0];
rz(-1.1739434) q[0];
sx q[0];
rz(-1.5482192) q[0];
rz(0.37823996) q[1];
sx q[1];
rz(-2.389237) q[1];
sx q[1];
rz(2.7517448) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5711576) q[0];
sx q[0];
rz(-1.3754554) q[0];
sx q[0];
rz(-0.51049149) q[0];
rz(-0.016640113) q[2];
sx q[2];
rz(-1.8748054) q[2];
sx q[2];
rz(1.2165516) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2521542) q[1];
sx q[1];
rz(-2.5233626) q[1];
sx q[1];
rz(-0.63207027) q[1];
rz(-pi) q[2];
rz(0.25501437) q[3];
sx q[3];
rz(-1.3552291) q[3];
sx q[3];
rz(0.6620342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0924015) q[2];
sx q[2];
rz(-1.1213877) q[2];
sx q[2];
rz(-1.8828877) q[2];
rz(-0.24426584) q[3];
sx q[3];
rz(-1.786307) q[3];
sx q[3];
rz(-0.99610966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8213537) q[0];
sx q[0];
rz(-1.9122253) q[0];
sx q[0];
rz(1.356333) q[0];
rz(2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(-0.43620268) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2552196) q[0];
sx q[0];
rz(-0.56111911) q[0];
sx q[0];
rz(-2.5291689) q[0];
rz(-pi) q[1];
rz(2.2460552) q[2];
sx q[2];
rz(-2.6288599) q[2];
sx q[2];
rz(0.89287478) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8456775) q[1];
sx q[1];
rz(-1.8009754) q[1];
sx q[1];
rz(-0.796411) q[1];
x q[2];
rz(-0.14472503) q[3];
sx q[3];
rz(-0.44071482) q[3];
sx q[3];
rz(-3.0851229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1205552) q[2];
sx q[2];
rz(-1.1641116) q[2];
sx q[2];
rz(0.54350054) q[2];
rz(3.0920658) q[3];
sx q[3];
rz(-1.4855569) q[3];
sx q[3];
rz(1.653695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5307584) q[0];
sx q[0];
rz(-0.13228358) q[0];
sx q[0];
rz(-0.079205967) q[0];
rz(2.7414956) q[1];
sx q[1];
rz(-2.2875417) q[1];
sx q[1];
rz(1.0265464) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75321046) q[0];
sx q[0];
rz(-2.3738656) q[0];
sx q[0];
rz(1.9110762) q[0];
rz(2.8978149) q[2];
sx q[2];
rz(-0.37636435) q[2];
sx q[2];
rz(1.7610904) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3803097) q[1];
sx q[1];
rz(-1.4766465) q[1];
sx q[1];
rz(0.32431079) q[1];
rz(1.6502569) q[3];
sx q[3];
rz(-1.9206408) q[3];
sx q[3];
rz(-0.81806483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93572271) q[2];
sx q[2];
rz(-1.8368072) q[2];
sx q[2];
rz(1.6892461) q[2];
rz(0.82990372) q[3];
sx q[3];
rz(-2.2645576) q[3];
sx q[3];
rz(2.1867866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5230539) q[0];
sx q[0];
rz(-1.1445615) q[0];
sx q[0];
rz(-1.6339697) q[0];
rz(-1.8745096) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(2.1562259) q[2];
sx q[2];
rz(-2.0337238) q[2];
sx q[2];
rz(0.82401333) q[2];
rz(2.6250056) q[3];
sx q[3];
rz(-2.5174375) q[3];
sx q[3];
rz(0.843912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
