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
rz(-0.77684075) q[0];
sx q[0];
rz(2.2651894) q[0];
sx q[0];
rz(9.6776008) q[0];
rz(-1.9934935) q[1];
sx q[1];
rz(-0.91528457) q[1];
sx q[1];
rz(-2.3372056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1348477) q[0];
sx q[0];
rz(-2.6383556) q[0];
sx q[0];
rz(0.87849599) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7115066) q[2];
sx q[2];
rz(-1.3580492) q[2];
sx q[2];
rz(1.2141922) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8411257) q[1];
sx q[1];
rz(-1.9411095) q[1];
sx q[1];
rz(1.200586) q[1];
rz(-pi) q[2];
rz(2.7976296) q[3];
sx q[3];
rz(-1.4620145) q[3];
sx q[3];
rz(0.0057084486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.09314166) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(0.87240458) q[2];
rz(2.8060272) q[3];
sx q[3];
rz(-1.395697) q[3];
sx q[3];
rz(0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5035079) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(0.22915325) q[0];
rz(-0.19042641) q[1];
sx q[1];
rz(-1.3399905) q[1];
sx q[1];
rz(-0.71358877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9035943) q[0];
sx q[0];
rz(-2.2269214) q[0];
sx q[0];
rz(-2.800992) q[0];
x q[1];
rz(2.3043556) q[2];
sx q[2];
rz(-1.1568489) q[2];
sx q[2];
rz(-0.6809823) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59871626) q[1];
sx q[1];
rz(-0.89756723) q[1];
sx q[1];
rz(-1.190541) q[1];
rz(-pi) q[2];
rz(-2.3768164) q[3];
sx q[3];
rz(-1.4524609) q[3];
sx q[3];
rz(1.8272994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4552292) q[2];
sx q[2];
rz(-2.8296622) q[2];
sx q[2];
rz(-0.4757821) q[2];
rz(2.0717715) q[3];
sx q[3];
rz(-2.1333623) q[3];
sx q[3];
rz(2.3750335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.0692724) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(-2.6909289) q[1];
sx q[1];
rz(-1.9602937) q[1];
sx q[1];
rz(2.252069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43972507) q[0];
sx q[0];
rz(-2.151818) q[0];
sx q[0];
rz(-1.0807316) q[0];
rz(-pi) q[1];
rz(2.6085147) q[2];
sx q[2];
rz(-2.0642082) q[2];
sx q[2];
rz(0.62668884) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3865859) q[1];
sx q[1];
rz(-2.4535547) q[1];
sx q[1];
rz(1.2351456) q[1];
rz(1.9185478) q[3];
sx q[3];
rz(-1.2158211) q[3];
sx q[3];
rz(2.6048911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.687261) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.8827776) q[2];
rz(1.9076294) q[3];
sx q[3];
rz(-1.1326658) q[3];
sx q[3];
rz(-3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4312129) q[0];
sx q[0];
rz(-0.90573913) q[0];
sx q[0];
rz(1.6718965) q[0];
rz(-1.9944893) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(-2.240644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9757864) q[0];
sx q[0];
rz(-2.8167509) q[0];
sx q[0];
rz(-0.27004098) q[0];
rz(-2.5261648) q[2];
sx q[2];
rz(-1.0503029) q[2];
sx q[2];
rz(-0.25582886) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.683953) q[1];
sx q[1];
rz(-1.6648653) q[1];
sx q[1];
rz(-0.72503083) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73293682) q[3];
sx q[3];
rz(-1.4519435) q[3];
sx q[3];
rz(2.262923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49742302) q[2];
sx q[2];
rz(-1.7449417) q[2];
sx q[2];
rz(-0.72745848) q[2];
rz(0.71508956) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(-2.6434744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46399507) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(-2.4053307) q[0];
rz(-2.5212506) q[1];
sx q[1];
rz(-2.5963929) q[1];
sx q[1];
rz(-1.5249407) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983561) q[0];
sx q[0];
rz(-1.5806206) q[0];
sx q[0];
rz(3.0863484) q[0];
x q[1];
rz(0.97648804) q[2];
sx q[2];
rz(-2.7433878) q[2];
sx q[2];
rz(-1.2334241) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.013629524) q[1];
sx q[1];
rz(-1.1971601) q[1];
sx q[1];
rz(0.84898265) q[1];
rz(-pi) q[2];
rz(-0.81186632) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(-0.52385274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6058495) q[2];
sx q[2];
rz(-2.3652786) q[2];
sx q[2];
rz(0.95174754) q[2];
rz(-2.7239299) q[3];
sx q[3];
rz(-2.8916736) q[3];
sx q[3];
rz(1.1550268) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8010537) q[0];
sx q[0];
rz(-1.8908353) q[0];
sx q[0];
rz(-0.75403768) q[0];
rz(-0.89318371) q[1];
sx q[1];
rz(-2.1008396) q[1];
sx q[1];
rz(2.975614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170068) q[0];
sx q[0];
rz(-2.4501194) q[0];
sx q[0];
rz(2.9777479) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50633035) q[2];
sx q[2];
rz(-1.4131318) q[2];
sx q[2];
rz(-2.8540996) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2365723) q[1];
sx q[1];
rz(-0.81099866) q[1];
sx q[1];
rz(0.6375957) q[1];
rz(-pi) q[2];
rz(-1.1096438) q[3];
sx q[3];
rz(-2.8361317) q[3];
sx q[3];
rz(2.4003254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37504998) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(-3.1362015) q[2];
rz(-0.088717669) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(2.3458792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2583112) q[0];
sx q[0];
rz(-1.2092051) q[0];
sx q[0];
rz(0.2051556) q[0];
rz(-0.908665) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(-0.044513449) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9712898) q[0];
sx q[0];
rz(-1.7006386) q[0];
sx q[0];
rz(1.4128039) q[0];
rz(-pi) q[1];
rz(-1.4304543) q[2];
sx q[2];
rz(-1.6082885) q[2];
sx q[2];
rz(0.692779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8537112) q[1];
sx q[1];
rz(-1.3559623) q[1];
sx q[1];
rz(-2.1532405) q[1];
rz(-pi) q[2];
rz(-1.7576587) q[3];
sx q[3];
rz(-3.0198041) q[3];
sx q[3];
rz(0.78021061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8160416) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(-2.6204056) q[2];
rz(-0.19736396) q[3];
sx q[3];
rz(-1.7557996) q[3];
sx q[3];
rz(0.50281966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407783) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(-0.25522301) q[0];
rz(-0.1420282) q[1];
sx q[1];
rz(-1.7250926) q[1];
sx q[1];
rz(2.7878888) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7434553) q[0];
sx q[0];
rz(-0.48851911) q[0];
sx q[0];
rz(-2.3228881) q[0];
x q[1];
rz(-2.2204705) q[2];
sx q[2];
rz(-2.4068458) q[2];
sx q[2];
rz(2.2237847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4434017) q[1];
sx q[1];
rz(-0.56487067) q[1];
sx q[1];
rz(2.2400212) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3925416) q[3];
sx q[3];
rz(-2.0045677) q[3];
sx q[3];
rz(1.449079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8346617) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(-1.6174512) q[2];
rz(-0.8664242) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(-2.5518937) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6397112) q[0];
sx q[0];
rz(-0.15637936) q[0];
sx q[0];
rz(1.3524652) q[0];
rz(1.9068708) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(-2.629705) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05737948) q[0];
sx q[0];
rz(-1.6239927) q[0];
sx q[0];
rz(1.4131111) q[0];
rz(-2.3761231) q[2];
sx q[2];
rz(-1.4569032) q[2];
sx q[2];
rz(-1.423045) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33495397) q[1];
sx q[1];
rz(-2.6158164) q[1];
sx q[1];
rz(-2.5277565) q[1];
rz(-pi) q[2];
rz(0.98661042) q[3];
sx q[3];
rz(-1.8711149) q[3];
sx q[3];
rz(2.9943313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0074244) q[2];
sx q[2];
rz(-1.3205426) q[2];
sx q[2];
rz(-2.7443938) q[2];
rz(-1.0154137) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(-2.5889682) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7427202) q[0];
sx q[0];
rz(-1.9341368) q[0];
sx q[0];
rz(-0.43750986) q[0];
rz(-2.4028026) q[1];
sx q[1];
rz(-0.80016017) q[1];
sx q[1];
rz(-1.4791666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9112944) q[0];
sx q[0];
rz(-2.0574967) q[0];
sx q[0];
rz(1.6337956) q[0];
x q[1];
rz(-0.11585856) q[2];
sx q[2];
rz(-1.2136974) q[2];
sx q[2];
rz(0.060572123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2228415) q[1];
sx q[1];
rz(-1.4645618) q[1];
sx q[1];
rz(1.5083471) q[1];
x q[2];
rz(0.71120925) q[3];
sx q[3];
rz(-2.2235302) q[3];
sx q[3];
rz(1.1389331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2207569) q[2];
sx q[2];
rz(-1.7831384) q[2];
sx q[2];
rz(0.1753359) q[2];
rz(1.0026503) q[3];
sx q[3];
rz(-2.9084539) q[3];
sx q[3];
rz(1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6709082) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(0.15057527) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(-0.19828037) q[2];
sx q[2];
rz(-1.5783327) q[2];
sx q[2];
rz(-0.21680149) q[2];
rz(2.3948432) q[3];
sx q[3];
rz(-2.8058553) q[3];
sx q[3];
rz(-0.6466858) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
