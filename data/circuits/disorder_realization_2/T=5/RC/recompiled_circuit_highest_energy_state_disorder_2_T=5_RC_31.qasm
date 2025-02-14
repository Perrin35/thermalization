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
rz(1.6588563) q[0];
sx q[0];
rz(-0.98134494) q[0];
sx q[0];
rz(-2.0438097) q[0];
rz(2.3535347) q[1];
sx q[1];
rz(-1.0507974) q[1];
sx q[1];
rz(0.0069590574) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46411447) q[0];
sx q[0];
rz(-0.65271806) q[0];
sx q[0];
rz(1.663289) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9240146) q[2];
sx q[2];
rz(-1.7597464) q[2];
sx q[2];
rz(1.762977) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3462249) q[1];
sx q[1];
rz(-0.2381033) q[1];
sx q[1];
rz(-1.7134929) q[1];
rz(-pi) q[2];
rz(2.801729) q[3];
sx q[3];
rz(-1.6714665) q[3];
sx q[3];
rz(-1.0235746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9030582) q[2];
sx q[2];
rz(-1.7944585) q[2];
sx q[2];
rz(-1.6713589) q[2];
rz(3.0155731) q[3];
sx q[3];
rz(-0.44243789) q[3];
sx q[3];
rz(-2.457705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0266492) q[0];
sx q[0];
rz(-2.6925955) q[0];
sx q[0];
rz(-2.5057416) q[0];
rz(2.3682829) q[1];
sx q[1];
rz(-0.67716235) q[1];
sx q[1];
rz(1.6962475) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4659368) q[0];
sx q[0];
rz(-2.0147062) q[0];
sx q[0];
rz(0.84966425) q[0];
rz(0.47451203) q[2];
sx q[2];
rz(-2.2998428) q[2];
sx q[2];
rz(-0.10708933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6856001) q[1];
sx q[1];
rz(-1.7280792) q[1];
sx q[1];
rz(0.70684163) q[1];
rz(-pi) q[2];
x q[2];
rz(0.068447114) q[3];
sx q[3];
rz(-1.7020149) q[3];
sx q[3];
rz(0.4438627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14629743) q[2];
sx q[2];
rz(-1.3714906) q[2];
sx q[2];
rz(-0.13776097) q[2];
rz(-0.45608258) q[3];
sx q[3];
rz(-2.5465953) q[3];
sx q[3];
rz(0.47541398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.544203) q[0];
sx q[0];
rz(-1.5578288) q[0];
sx q[0];
rz(0.57918817) q[0];
rz(0.10313615) q[1];
sx q[1];
rz(-1.3214279) q[1];
sx q[1];
rz(0.78027049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1121169) q[0];
sx q[0];
rz(-1.5896086) q[0];
sx q[0];
rz(0.095186724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4109072) q[2];
sx q[2];
rz(-0.89582755) q[2];
sx q[2];
rz(-2.8596148) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0611813) q[1];
sx q[1];
rz(-1.1763417) q[1];
sx q[1];
rz(-2.9030672) q[1];
rz(-pi) q[2];
rz(-1.8336201) q[3];
sx q[3];
rz(-0.71280957) q[3];
sx q[3];
rz(-1.7607016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66037336) q[2];
sx q[2];
rz(-1.6796835) q[2];
sx q[2];
rz(-3.050728) q[2];
rz(1.4800492) q[3];
sx q[3];
rz(-1.816498) q[3];
sx q[3];
rz(2.3265694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.018983) q[0];
sx q[0];
rz(-0.18773395) q[0];
sx q[0];
rz(-2.6222099) q[0];
rz(1.9620365) q[1];
sx q[1];
rz(-2.4603381) q[1];
sx q[1];
rz(0.048351668) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40019401) q[0];
sx q[0];
rz(-2.1693008) q[0];
sx q[0];
rz(-2.5116176) q[0];
rz(-pi) q[1];
rz(-2.2273793) q[2];
sx q[2];
rz(-0.58744741) q[2];
sx q[2];
rz(2.812831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4614233) q[1];
sx q[1];
rz(-1.5986414) q[1];
sx q[1];
rz(2.6538886) q[1];
rz(-pi) q[2];
rz(-2.8818733) q[3];
sx q[3];
rz(-0.81213299) q[3];
sx q[3];
rz(-0.5838697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68253303) q[2];
sx q[2];
rz(-2.8352663) q[2];
sx q[2];
rz(1.4502067) q[2];
rz(-2.0356483) q[3];
sx q[3];
rz(-1.9927497) q[3];
sx q[3];
rz(-0.86627427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3806216) q[0];
sx q[0];
rz(-0.80046099) q[0];
sx q[0];
rz(-2.3642484) q[0];
rz(2.6457973) q[1];
sx q[1];
rz(-0.40758857) q[1];
sx q[1];
rz(2.9487603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042943311) q[0];
sx q[0];
rz(-1.0355562) q[0];
sx q[0];
rz(-1.213206) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76754153) q[2];
sx q[2];
rz(-1.3965675) q[2];
sx q[2];
rz(-1.4105287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1950394) q[1];
sx q[1];
rz(-1.6583539) q[1];
sx q[1];
rz(1.8680496) q[1];
rz(-pi) q[2];
rz(-2.9787872) q[3];
sx q[3];
rz(-1.9678402) q[3];
sx q[3];
rz(0.59418488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.58482802) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(-0.83621109) q[2];
rz(-2.1861475) q[3];
sx q[3];
rz(-1.4232057) q[3];
sx q[3];
rz(0.53171617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504976) q[0];
sx q[0];
rz(-1.9170772) q[0];
sx q[0];
rz(-0.033705458) q[0];
rz(0.91824245) q[1];
sx q[1];
rz(-2.4372209) q[1];
sx q[1];
rz(-0.69923002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3501773) q[0];
sx q[0];
rz(-0.79424131) q[0];
sx q[0];
rz(2.3167531) q[0];
x q[1];
rz(1.6246454) q[2];
sx q[2];
rz(-1.3837746) q[2];
sx q[2];
rz(2.5973926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67656006) q[1];
sx q[1];
rz(-2.0102215) q[1];
sx q[1];
rz(-1.7656754) q[1];
rz(-2.8630303) q[3];
sx q[3];
rz(-1.2435561) q[3];
sx q[3];
rz(2.183941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0596727) q[2];
sx q[2];
rz(-2.4884188) q[2];
sx q[2];
rz(-3.0461404) q[2];
rz(-2.0898315) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(2.2473647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8555701) q[0];
sx q[0];
rz(-0.48208553) q[0];
sx q[0];
rz(-1.9739738) q[0];
rz(2.7883912) q[1];
sx q[1];
rz(-0.51274061) q[1];
sx q[1];
rz(-0.28164992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324556) q[0];
sx q[0];
rz(-1.0971709) q[0];
sx q[0];
rz(3.043463) q[0];
rz(1.4640305) q[2];
sx q[2];
rz(-2.2153478) q[2];
sx q[2];
rz(2.0951955) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6503164) q[1];
sx q[1];
rz(-1.4230844) q[1];
sx q[1];
rz(3.0056535) q[1];
x q[2];
rz(0.8696712) q[3];
sx q[3];
rz(-0.68415776) q[3];
sx q[3];
rz(1.9811833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99792751) q[2];
sx q[2];
rz(-1.2219656) q[2];
sx q[2];
rz(-1.3227051) q[2];
rz(-1.6392684) q[3];
sx q[3];
rz(-1.084525) q[3];
sx q[3];
rz(3.0433906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1439576) q[0];
sx q[0];
rz(-1.8956381) q[0];
sx q[0];
rz(-2.8676046) q[0];
rz(-0.59448376) q[1];
sx q[1];
rz(-1.4130054) q[1];
sx q[1];
rz(-0.53057539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0230867) q[0];
sx q[0];
rz(-0.039561633) q[0];
sx q[0];
rz(0.096574144) q[0];
rz(3.0184348) q[2];
sx q[2];
rz(-0.55610699) q[2];
sx q[2];
rz(-1.7983939) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6484687) q[1];
sx q[1];
rz(-1.1965355) q[1];
sx q[1];
rz(1.7497653) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.382393) q[3];
sx q[3];
rz(-0.72160463) q[3];
sx q[3];
rz(2.5327275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8990367) q[2];
sx q[2];
rz(-1.2071004) q[2];
sx q[2];
rz(0.25699082) q[2];
rz(-1.9422003) q[3];
sx q[3];
rz(-1.896984) q[3];
sx q[3];
rz(-1.6283584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5892107) q[0];
sx q[0];
rz(-0.95443812) q[0];
sx q[0];
rz(-2.6115665) q[0];
rz(1.5026622) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(-0.93856215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6185222) q[0];
sx q[0];
rz(-1.1022727) q[0];
sx q[0];
rz(2.4836088) q[0];
rz(-pi) q[1];
x q[1];
rz(1.45088) q[2];
sx q[2];
rz(-1.7133822) q[2];
sx q[2];
rz(0.8220807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1129866) q[1];
sx q[1];
rz(-0.48278207) q[1];
sx q[1];
rz(-1.0989936) q[1];
rz(-pi) q[2];
rz(-3.0051232) q[3];
sx q[3];
rz(-2.1519063) q[3];
sx q[3];
rz(-2.7032397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75371257) q[2];
sx q[2];
rz(-2.2215999) q[2];
sx q[2];
rz(-2.06125) q[2];
rz(-2.3267817) q[3];
sx q[3];
rz(-1.9459414) q[3];
sx q[3];
rz(-1.7241534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.998488) q[0];
sx q[0];
rz(-1.4838706) q[0];
sx q[0];
rz(2.0527573) q[0];
rz(3.0740652) q[1];
sx q[1];
rz(-1.8636401) q[1];
sx q[1];
rz(-0.75928226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5048333) q[0];
sx q[0];
rz(-0.28226109) q[0];
sx q[0];
rz(1.1310079) q[0];
rz(-pi) q[1];
rz(-1.4467054) q[2];
sx q[2];
rz(-1.795009) q[2];
sx q[2];
rz(-1.8607163) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2226505) q[1];
sx q[1];
rz(-0.83178751) q[1];
sx q[1];
rz(1.8636501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5645026) q[3];
sx q[3];
rz(-2.7365757) q[3];
sx q[3];
rz(-1.7187723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33409432) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(2.5661772) q[2];
rz(2.5790162) q[3];
sx q[3];
rz(-0.96499363) q[3];
sx q[3];
rz(-1.3728728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0811049) q[0];
sx q[0];
rz(-1.8555547) q[0];
sx q[0];
rz(-0.9077358) q[0];
rz(-2.0545215) q[1];
sx q[1];
rz(-1.1320976) q[1];
sx q[1];
rz(-0.017398106) q[1];
rz(1.4292064) q[2];
sx q[2];
rz(-1.7238486) q[2];
sx q[2];
rz(0.42949745) q[2];
rz(-1.1997052) q[3];
sx q[3];
rz(-1.0047079) q[3];
sx q[3];
rz(-2.9550002) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
