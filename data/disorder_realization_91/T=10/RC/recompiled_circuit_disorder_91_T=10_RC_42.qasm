OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(-0.17106549) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(3.4847335) q[1];
sx q[1];
rz(7.6141678) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5509697) q[0];
sx q[0];
rz(-1.5184214) q[0];
sx q[0];
rz(2.4916324) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6257329) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(3.0992103) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0099437873) q[1];
sx q[1];
rz(-2.8899) q[1];
sx q[1];
rz(-1.6394872) q[1];
x q[2];
rz(2.8475548) q[3];
sx q[3];
rz(-2.4774744) q[3];
sx q[3];
rz(2.3703863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3657637) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(1.8120871) q[2];
rz(-1.7261516) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(1.6888899) q[0];
rz(0.20092043) q[1];
sx q[1];
rz(-1.0849489) q[1];
sx q[1];
rz(-0.30028775) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0694885) q[0];
sx q[0];
rz(-0.61394962) q[0];
sx q[0];
rz(-2.9314562) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47313182) q[2];
sx q[2];
rz(-2.5169249) q[2];
sx q[2];
rz(-1.9249137) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2806432) q[1];
sx q[1];
rz(-1.7935392) q[1];
sx q[1];
rz(1.2692979) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85864752) q[3];
sx q[3];
rz(-2.894069) q[3];
sx q[3];
rz(1.1502707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.369027) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(2.1035813) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5791941) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(0.38811362) q[0];
rz(-0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-0.31633502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4916723) q[0];
sx q[0];
rz(-2.8965817) q[0];
sx q[0];
rz(1.578376) q[0];
x q[1];
rz(2.5306273) q[2];
sx q[2];
rz(-0.66662153) q[2];
sx q[2];
rz(1.9463469) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9005147) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(-0.48836744) q[1];
rz(0.86795904) q[3];
sx q[3];
rz(-0.43324019) q[3];
sx q[3];
rz(0.058335282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1091653) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(2.9807828) q[2];
rz(-3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2407103) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(1.6279189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8320223) q[0];
sx q[0];
rz(-2.4674468) q[0];
sx q[0];
rz(3.0679697) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.081110031) q[2];
sx q[2];
rz(-2.1334279) q[2];
sx q[2];
rz(-0.62603355) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1308243) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(0.018647714) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4098445) q[3];
sx q[3];
rz(-1.2536612) q[3];
sx q[3];
rz(-1.5750615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.224219) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(-0.52880374) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(1.303724) q[0];
rz(2.6858792) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(-0.051503332) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951552) q[0];
sx q[0];
rz(-1.6704541) q[0];
sx q[0];
rz(-1.6874203) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0837469) q[2];
sx q[2];
rz(-2.8205928) q[2];
sx q[2];
rz(0.1333065) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0609378) q[1];
sx q[1];
rz(-2.6662711) q[1];
sx q[1];
rz(1.2259543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4476942) q[3];
sx q[3];
rz(-2.3007292) q[3];
sx q[3];
rz(2.6134932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6872528) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(-2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0705868) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(0.083267033) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986388) q[0];
sx q[0];
rz(-2.8790701) q[0];
sx q[0];
rz(-0.72857626) q[0];
rz(-pi) q[1];
rz(1.8385356) q[2];
sx q[2];
rz(-1.3053075) q[2];
sx q[2];
rz(0.43648411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6492918) q[1];
sx q[1];
rz(-1.2052844) q[1];
sx q[1];
rz(3.027012) q[1];
rz(-1.2721328) q[3];
sx q[3];
rz(-2.1080058) q[3];
sx q[3];
rz(0.99029535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0319556) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(2.8921195) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(1.680254) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-2.0239963) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5516978) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(-2.9445573) q[0];
rz(-pi) q[1];
rz(-0.20861161) q[2];
sx q[2];
rz(-0.87000712) q[2];
sx q[2];
rz(-2.029062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1114137) q[1];
sx q[1];
rz(-1.3707268) q[1];
sx q[1];
rz(2.5740037) q[1];
rz(-pi) q[2];
rz(-0.53447978) q[3];
sx q[3];
rz(-1.3325053) q[3];
sx q[3];
rz(-2.0508545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(2.7548742) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(-2.1462671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052208386) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(2.9242933) q[0];
rz(3.1106588) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-2.4826179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17651672) q[0];
sx q[0];
rz(-2.342431) q[0];
sx q[0];
rz(1.6060711) q[0];
rz(-0.89491567) q[2];
sx q[2];
rz(-2.1721828) q[2];
sx q[2];
rz(-0.88047699) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3498889) q[1];
sx q[1];
rz(-0.83194299) q[1];
sx q[1];
rz(1.2757343) q[1];
rz(-2.213845) q[3];
sx q[3];
rz(-1.6258996) q[3];
sx q[3];
rz(0.4004713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7028246) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(0.32021114) q[2];
rz(-1.6163588) q[3];
sx q[3];
rz(-1.1956918) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76072389) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(-1.0702417) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.9314996) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4142128) q[0];
sx q[0];
rz(-2.7276037) q[0];
sx q[0];
rz(2.1075641) q[0];
rz(-0.25173431) q[2];
sx q[2];
rz(-0.98099698) q[2];
sx q[2];
rz(0.41433197) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.853211) q[1];
sx q[1];
rz(-2.2182811) q[1];
sx q[1];
rz(-0.37154571) q[1];
x q[2];
rz(-0.61981598) q[3];
sx q[3];
rz(-1.8780939) q[3];
sx q[3];
rz(0.67051552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(-0.56662095) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728445) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.7096827) q[0];
rz(0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-0.79968232) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0737338) q[0];
sx q[0];
rz(-1.6556544) q[0];
sx q[0];
rz(1.2220864) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4893555) q[2];
sx q[2];
rz(-1.4105721) q[2];
sx q[2];
rz(-0.83934957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77196808) q[1];
sx q[1];
rz(-0.63056417) q[1];
sx q[1];
rz(-2.1715013) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47761376) q[3];
sx q[3];
rz(-0.16371809) q[3];
sx q[3];
rz(-1.175566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(0.69520673) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
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
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.862539) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(-2.8056801) q[2];
sx q[2];
rz(-1.5699785) q[2];
sx q[2];
rz(1.7287398) q[2];
rz(-2.3435076) q[3];
sx q[3];
rz(-1.4555664) q[3];
sx q[3];
rz(1.2496787) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];