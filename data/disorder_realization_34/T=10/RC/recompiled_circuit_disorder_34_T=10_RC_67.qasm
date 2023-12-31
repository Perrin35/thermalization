OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(5.8435506) q[0];
sx q[0];
rz(6.2018659) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5188132) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(1.326484) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15687234) q[2];
sx q[2];
rz(-1.9777159) q[2];
sx q[2];
rz(-3.0770609) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70071917) q[1];
sx q[1];
rz(-2.6881725) q[1];
sx q[1];
rz(-1.8666301) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8257915) q[3];
sx q[3];
rz(-1.485802) q[3];
sx q[3];
rz(2.1080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2471182) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(-1.5995021) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2540934) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(0.37503606) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(2.9017752) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46562425) q[0];
sx q[0];
rz(-1.3249319) q[0];
sx q[0];
rz(-1.7322391) q[0];
x q[1];
rz(-1.2698783) q[2];
sx q[2];
rz(-1.3687203) q[2];
sx q[2];
rz(-2.2965455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5259243) q[1];
sx q[1];
rz(-0.79155542) q[1];
sx q[1];
rz(0.073992373) q[1];
rz(-0.85571839) q[3];
sx q[3];
rz(-0.021589605) q[3];
sx q[3];
rz(1.74687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(-1.3416393) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-2.4734316) q[0];
rz(-1.4913303) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6324368) q[0];
sx q[0];
rz(-1.3226489) q[0];
sx q[0];
rz(-3.0993673) q[0];
rz(2.9124444) q[2];
sx q[2];
rz(-2.5742968) q[2];
sx q[2];
rz(-1.0457525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26820688) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(0.98209776) q[1];
rz(-pi) q[2];
rz(3.021574) q[3];
sx q[3];
rz(-2.2285322) q[3];
sx q[3];
rz(-2.0002055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(0.032547396) q[2];
rz(-0.35999808) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(-1.6548086) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2457755) q[0];
sx q[0];
rz(-2.7005807) q[0];
sx q[0];
rz(-2.4367024) q[0];
rz(-2.5321132) q[2];
sx q[2];
rz(-1.8595427) q[2];
sx q[2];
rz(0.8537054) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6390037) q[1];
sx q[1];
rz(-1.2856312) q[1];
sx q[1];
rz(1.7267978) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5634469) q[3];
sx q[3];
rz(-2.3971933) q[3];
sx q[3];
rz(-3.0878382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(-2.0452943) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(-1.6812356) q[0];
rz(-1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(-0.016074093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0468633) q[0];
sx q[0];
rz(-1.0445347) q[0];
sx q[0];
rz(2.2577283) q[0];
rz(-pi) q[1];
rz(1.5526505) q[2];
sx q[2];
rz(-2.7108253) q[2];
sx q[2];
rz(-0.22242966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4453455) q[1];
sx q[1];
rz(-1.0187341) q[1];
sx q[1];
rz(0.50260431) q[1];
rz(-pi) q[2];
rz(3.1297242) q[3];
sx q[3];
rz(-2.7792645) q[3];
sx q[3];
rz(-3.0167937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(-1.0970998) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(-2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(0.81470195) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(1.9168568) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46366102) q[0];
sx q[0];
rz(-1.3555129) q[0];
sx q[0];
rz(-0.93925516) q[0];
rz(-pi) q[1];
rz(0.20829006) q[2];
sx q[2];
rz(-2.0049094) q[2];
sx q[2];
rz(1.2693894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9310301) q[1];
sx q[1];
rz(-1.5157053) q[1];
sx q[1];
rz(-0.964826) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8858658) q[3];
sx q[3];
rz(-0.41318196) q[3];
sx q[3];
rz(-0.43858389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(0.93969807) q[2];
rz(0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796824) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(0.58037037) q[0];
rz(1.0549818) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(-2.4408128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66877767) q[0];
sx q[0];
rz(-2.6758709) q[0];
sx q[0];
rz(-1.3186243) q[0];
rz(2.8009311) q[2];
sx q[2];
rz(-0.33764687) q[2];
sx q[2];
rz(1.1866736) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5138445) q[1];
sx q[1];
rz(-2.5690837) q[1];
sx q[1];
rz(-1.7379012) q[1];
x q[2];
rz(-0.084214597) q[3];
sx q[3];
rz(-2.299752) q[3];
sx q[3];
rz(0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-3.11943) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.78684029) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4448924) q[0];
sx q[0];
rz(-0.80388821) q[0];
sx q[0];
rz(-0.21960396) q[0];
rz(-pi) q[1];
rz(-2.5052091) q[2];
sx q[2];
rz(-2.4578265) q[2];
sx q[2];
rz(-2.6725339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29854952) q[1];
sx q[1];
rz(-1.9798568) q[1];
sx q[1];
rz(-2.4119853) q[1];
rz(-2.2303921) q[3];
sx q[3];
rz(-1.824114) q[3];
sx q[3];
rz(2.5384056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(0.78424224) q[2];
rz(0.50576058) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(-2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.4760251) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(-1.7766215) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7576335) q[0];
sx q[0];
rz(-0.74476349) q[0];
sx q[0];
rz(-0.061116771) q[0];
x q[1];
rz(1.1599837) q[2];
sx q[2];
rz(-1.4940726) q[2];
sx q[2];
rz(0.086417925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0627928) q[1];
sx q[1];
rz(-1.2275057) q[1];
sx q[1];
rz(2.417056) q[1];
rz(-pi) q[2];
rz(-1.8678719) q[3];
sx q[3];
rz(-1.3798514) q[3];
sx q[3];
rz(-1.8553268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1086796) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(0.03406295) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0868527) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(1.6636794) q[0];
rz(2.058303) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(2.1733984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6598845) q[0];
sx q[0];
rz(-1.1904926) q[0];
sx q[0];
rz(-0.36614059) q[0];
rz(-0.88629006) q[2];
sx q[2];
rz(-0.61294014) q[2];
sx q[2];
rz(0.71000368) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89322461) q[1];
sx q[1];
rz(-1.4239422) q[1];
sx q[1];
rz(2.8180772) q[1];
rz(-0.00023437436) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(-2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0782464) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(-0.001312288) q[2];
rz(1.1408268) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.39682) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(0.36623476) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(-1.2583854) q[2];
sx q[2];
rz(-0.84597833) q[2];
sx q[2];
rz(-0.99920263) q[2];
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
