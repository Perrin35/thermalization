OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1333753) q[0];
sx q[0];
rz(-1.8741338) q[0];
sx q[0];
rz(-3.128669) q[0];
rz(0.68459964) q[1];
sx q[1];
rz(-2.3426988) q[1];
sx q[1];
rz(1.0577143) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2051281) q[0];
sx q[0];
rz(-1.797344) q[0];
sx q[0];
rz(-0.675987) q[0];
x q[1];
rz(1.0619668) q[2];
sx q[2];
rz(-2.2812009) q[2];
sx q[2];
rz(1.8343385) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.31725059) q[1];
sx q[1];
rz(-1.2095272) q[1];
sx q[1];
rz(0.74202219) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1790479) q[3];
sx q[3];
rz(-1.5314252) q[3];
sx q[3];
rz(-1.5567428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58618033) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(0.093322873) q[2];
rz(-3.1210476) q[3];
sx q[3];
rz(-2.8829657) q[3];
sx q[3];
rz(-1.3625905) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697295) q[0];
sx q[0];
rz(-1.7427895) q[0];
sx q[0];
rz(-0.82759696) q[0];
rz(2.1780275) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(2.4049984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4086491) q[0];
sx q[0];
rz(-1.8445065) q[0];
sx q[0];
rz(3.0628171) q[0];
x q[1];
rz(-2.170855) q[2];
sx q[2];
rz(-1.0781556) q[2];
sx q[2];
rz(-0.095183177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.40562661) q[1];
sx q[1];
rz(-2.7752004) q[1];
sx q[1];
rz(-0.97586164) q[1];
rz(0.16421825) q[3];
sx q[3];
rz(-3.0053557) q[3];
sx q[3];
rz(0.0020023684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5654512) q[2];
sx q[2];
rz(-2.5665923) q[2];
sx q[2];
rz(-2.3973993) q[2];
rz(0.60892504) q[3];
sx q[3];
rz(-2.3600793) q[3];
sx q[3];
rz(-1.6412546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.7805444) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(1.3737099) q[0];
rz(0.79958493) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(1.132157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4904259) q[0];
sx q[0];
rz(-1.4803783) q[0];
sx q[0];
rz(-1.3549442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5930873) q[2];
sx q[2];
rz(-1.4906851) q[2];
sx q[2];
rz(-1.3673283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.76946041) q[1];
sx q[1];
rz(-2.4677708) q[1];
sx q[1];
rz(0.12427434) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6145094) q[3];
sx q[3];
rz(-0.96230405) q[3];
sx q[3];
rz(0.19245806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75480294) q[2];
sx q[2];
rz(-0.58480442) q[2];
sx q[2];
rz(1.7542138) q[2];
rz(2.7925708) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(2.6121228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.3997407) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(-0.35811785) q[0];
rz(2.070836) q[1];
sx q[1];
rz(-2.4929969) q[1];
sx q[1];
rz(1.4395641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0107797) q[0];
sx q[0];
rz(-2.3284916) q[0];
sx q[0];
rz(1.9365501) q[0];
x q[1];
rz(-2.1522572) q[2];
sx q[2];
rz(-0.78063595) q[2];
sx q[2];
rz(-2.6065741) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6144952) q[1];
sx q[1];
rz(-1.004809) q[1];
sx q[1];
rz(2.4742592) q[1];
rz(0.26243383) q[3];
sx q[3];
rz(-2.35733) q[3];
sx q[3];
rz(1.5465496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4293356) q[2];
sx q[2];
rz(-1.9108994) q[2];
sx q[2];
rz(-2.5168391) q[2];
rz(-1.6338927) q[3];
sx q[3];
rz(-0.75993901) q[3];
sx q[3];
rz(2.0190575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7655012) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(1.1464024) q[0];
rz(1.435185) q[1];
sx q[1];
rz(-2.621666) q[1];
sx q[1];
rz(0.82040876) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489481) q[0];
sx q[0];
rz(-1.2995509) q[0];
sx q[0];
rz(3.0503037) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90419714) q[2];
sx q[2];
rz(-1.7115286) q[2];
sx q[2];
rz(-2.4183969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.31215224) q[1];
sx q[1];
rz(-1.3669126) q[1];
sx q[1];
rz(-1.7704493) q[1];
rz(-pi) q[2];
rz(3.0919364) q[3];
sx q[3];
rz(-0.59030246) q[3];
sx q[3];
rz(-1.717339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4121805) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(-0.5385651) q[2];
rz(1.6719079) q[3];
sx q[3];
rz(-1.6939751) q[3];
sx q[3];
rz(-0.058549747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84022123) q[0];
sx q[0];
rz(-0.138962) q[0];
sx q[0];
rz(3.050991) q[0];
rz(0.36965707) q[1];
sx q[1];
rz(-1.5934817) q[1];
sx q[1];
rz(-0.37818092) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7525258) q[0];
sx q[0];
rz(-1.7992668) q[0];
sx q[0];
rz(0.48918251) q[0];
x q[1];
rz(1.6126851) q[2];
sx q[2];
rz(-1.48078) q[2];
sx q[2];
rz(2.5652094) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0269811) q[1];
sx q[1];
rz(-1.0798732) q[1];
sx q[1];
rz(1.4308962) q[1];
rz(-1.380081) q[3];
sx q[3];
rz(-1.6487953) q[3];
sx q[3];
rz(-1.5275148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6211264) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(0.10144083) q[2];
rz(1.2927879) q[3];
sx q[3];
rz(-2.3088876) q[3];
sx q[3];
rz(-1.4147991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
rz(1.3066278) q[0];
sx q[0];
rz(-1.2766301) q[0];
sx q[0];
rz(0.90245885) q[0];
rz(2.9023671) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(0.36144027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1162565) q[0];
sx q[0];
rz(-2.6317208) q[0];
sx q[0];
rz(-1.8007397) q[0];
x q[1];
rz(-2.7290384) q[2];
sx q[2];
rz(-0.13880402) q[2];
sx q[2];
rz(-2.1766162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4415008) q[1];
sx q[1];
rz(-0.48466408) q[1];
sx q[1];
rz(1.086471) q[1];
rz(2.802235) q[3];
sx q[3];
rz(-1.4213143) q[3];
sx q[3];
rz(0.71163346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0685737) q[2];
sx q[2];
rz(-1.34015) q[2];
sx q[2];
rz(2.2646591) q[2];
rz(2.1584611) q[3];
sx q[3];
rz(-1.5638331) q[3];
sx q[3];
rz(0.94299281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32901949) q[0];
sx q[0];
rz(-2.946377) q[0];
sx q[0];
rz(0.83874291) q[0];
rz(0.70873952) q[1];
sx q[1];
rz(-0.63116169) q[1];
sx q[1];
rz(-0.90604025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86676407) q[0];
sx q[0];
rz(-1.5497009) q[0];
sx q[0];
rz(0.72480185) q[0];
rz(-pi) q[1];
rz(-1.0308517) q[2];
sx q[2];
rz(-2.738852) q[2];
sx q[2];
rz(-0.38322057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38323572) q[1];
sx q[1];
rz(-0.89665186) q[1];
sx q[1];
rz(-0.70887776) q[1];
x q[2];
rz(1.2165478) q[3];
sx q[3];
rz(-0.63874001) q[3];
sx q[3];
rz(0.50001345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9771542) q[2];
sx q[2];
rz(-0.70455569) q[2];
sx q[2];
rz(2.0257115) q[2];
rz(-2.5130533) q[3];
sx q[3];
rz(-1.081859) q[3];
sx q[3];
rz(-0.61454296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3299265) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(2.9803168) q[0];
rz(-2.2992086) q[1];
sx q[1];
rz(-1.7120275) q[1];
sx q[1];
rz(2.9715723) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1415213) q[0];
sx q[0];
rz(-1.562644) q[0];
sx q[0];
rz(-3.1233112) q[0];
rz(-pi) q[1];
rz(-0.44012897) q[2];
sx q[2];
rz(-1.7227731) q[2];
sx q[2];
rz(0.51896799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5698118) q[1];
sx q[1];
rz(-1.7649635) q[1];
sx q[1];
rz(-1.1295736) q[1];
rz(-1.6415855) q[3];
sx q[3];
rz(-2.7578691) q[3];
sx q[3];
rz(0.078670382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9876447) q[2];
sx q[2];
rz(-1.5445856) q[2];
sx q[2];
rz(-1.6551931) q[2];
rz(-3.0814643) q[3];
sx q[3];
rz(-2.3190053) q[3];
sx q[3];
rz(-1.7000343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104367) q[0];
sx q[0];
rz(-0.031351723) q[0];
sx q[0];
rz(1.4309058) q[0];
rz(0.36262861) q[1];
sx q[1];
rz(-1.6094094) q[1];
sx q[1];
rz(1.3154202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5820438) q[0];
sx q[0];
rz(-1.8616195) q[0];
sx q[0];
rz(-1.8946339) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1926629) q[2];
sx q[2];
rz(-1.6844498) q[2];
sx q[2];
rz(-3.0494351) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9818287) q[1];
sx q[1];
rz(-1.8366645) q[1];
sx q[1];
rz(-0.60536107) q[1];
rz(-pi) q[2];
rz(-2.0995448) q[3];
sx q[3];
rz(-2.3662851) q[3];
sx q[3];
rz(0.14107832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0843087) q[2];
sx q[2];
rz(-1.5466651) q[2];
sx q[2];
rz(-0.16051897) q[2];
rz(2.6594243) q[3];
sx q[3];
rz(-0.67626685) q[3];
sx q[3];
rz(0.67960656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01874825) q[0];
sx q[0];
rz(-1.7457122) q[0];
sx q[0];
rz(2.2298298) q[0];
rz(-1.1625166) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(-0.087333655) q[2];
sx q[2];
rz(-1.4764016) q[2];
sx q[2];
rz(1.7511677) q[2];
rz(0.58070498) q[3];
sx q[3];
rz(-2.28149) q[3];
sx q[3];
rz(-2.1547439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
