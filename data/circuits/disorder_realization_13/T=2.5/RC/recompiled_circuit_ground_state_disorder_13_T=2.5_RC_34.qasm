OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2406727) q[0];
sx q[0];
rz(-2.9958041) q[0];
sx q[0];
rz(-1.2989651) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(-1.8008404) q[1];
sx q[1];
rz(0.10032108) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.318905) q[0];
sx q[0];
rz(-1.6226557) q[0];
sx q[0];
rz(2.479631) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47122987) q[2];
sx q[2];
rz(-1.8719851) q[2];
sx q[2];
rz(2.9848841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4397706) q[1];
sx q[1];
rz(-1.7793097) q[1];
sx q[1];
rz(-0.28065248) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2572977) q[3];
sx q[3];
rz(-2.1380205) q[3];
sx q[3];
rz(-0.49916609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4002865) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(-0.8068223) q[2];
rz(2.412879) q[3];
sx q[3];
rz(-2.3829134) q[3];
sx q[3];
rz(-1.2046643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.62419409) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(0.30701315) q[0];
rz(1.2767876) q[1];
sx q[1];
rz(-2.0025608) q[1];
sx q[1];
rz(2.0397287) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1043423) q[0];
sx q[0];
rz(-2.5581963) q[0];
sx q[0];
rz(-1.2520381) q[0];
x q[1];
rz(1.6346143) q[2];
sx q[2];
rz(-1.391727) q[2];
sx q[2];
rz(-0.049953559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4408177) q[1];
sx q[1];
rz(-1.3937181) q[1];
sx q[1];
rz(3.0040279) q[1];
x q[2];
rz(2.4073007) q[3];
sx q[3];
rz(-2.2059388) q[3];
sx q[3];
rz(1.7226302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.028194204) q[2];
sx q[2];
rz(-2.5598309) q[2];
sx q[2];
rz(0.096435189) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(-2.4436387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518799) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(-0.40019792) q[0];
rz(1.5125037) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(2.8834744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54329693) q[0];
sx q[0];
rz(-0.27320293) q[0];
sx q[0];
rz(2.1092723) q[0];
rz(1.6506972) q[2];
sx q[2];
rz(-1.6703484) q[2];
sx q[2];
rz(-2.4308506) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26545721) q[1];
sx q[1];
rz(-1.5837513) q[1];
sx q[1];
rz(1.575768) q[1];
x q[2];
rz(-2.5260365) q[3];
sx q[3];
rz(-2.8488692) q[3];
sx q[3];
rz(0.20197257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.021598024) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(0.032133948) q[2];
rz(-0.26767996) q[3];
sx q[3];
rz(-1.6240424) q[3];
sx q[3];
rz(-1.5816429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716008) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(0.084240325) q[0];
rz(-3.1056504) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(2.8004004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223136) q[0];
sx q[0];
rz(-1.8321165) q[0];
sx q[0];
rz(2.4909291) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58653763) q[2];
sx q[2];
rz(-1.3549202) q[2];
sx q[2];
rz(-0.12556533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2349639) q[1];
sx q[1];
rz(-0.90032691) q[1];
sx q[1];
rz(2.5690298) q[1];
rz(-0.59935948) q[3];
sx q[3];
rz(-1.1903569) q[3];
sx q[3];
rz(-0.75293702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42698947) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(-1.535447) q[2];
rz(-2.8793907) q[3];
sx q[3];
rz(-1.6180792) q[3];
sx q[3];
rz(-2.1867627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3839805) q[0];
sx q[0];
rz(-2.7181427) q[0];
sx q[0];
rz(-2.0641932) q[0];
rz(-2.7491838) q[1];
sx q[1];
rz(-3.0632186) q[1];
sx q[1];
rz(-2.1108625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01973132) q[0];
sx q[0];
rz(-1.1543589) q[0];
sx q[0];
rz(-2.5520578) q[0];
rz(-0.18056492) q[2];
sx q[2];
rz(-2.0624277) q[2];
sx q[2];
rz(0.59300834) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.64068078) q[1];
sx q[1];
rz(-1.7927153) q[1];
sx q[1];
rz(-2.9725916) q[1];
rz(-pi) q[2];
rz(-0.34319539) q[3];
sx q[3];
rz(-1.0846018) q[3];
sx q[3];
rz(-0.23517683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3173759) q[2];
sx q[2];
rz(-2.4874096) q[2];
sx q[2];
rz(0.83924323) q[2];
rz(0.9134891) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(-0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4588673) q[0];
sx q[0];
rz(-2.904628) q[0];
sx q[0];
rz(1.4759901) q[0];
rz(-2.7492375) q[1];
sx q[1];
rz(-1.0959492) q[1];
sx q[1];
rz(-2.5501693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8845733) q[0];
sx q[0];
rz(-1.9396175) q[0];
sx q[0];
rz(0.9414282) q[0];
x q[1];
rz(0.5742091) q[2];
sx q[2];
rz(-2.817135) q[2];
sx q[2];
rz(2.7221219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0589231) q[1];
sx q[1];
rz(-1.6367404) q[1];
sx q[1];
rz(-2.0772499) q[1];
rz(-pi) q[2];
rz(1.1272085) q[3];
sx q[3];
rz(-0.74723703) q[3];
sx q[3];
rz(-2.0153449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.300294) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(0.57233125) q[2];
rz(-2.9193997) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(-2.4967994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7543024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(-0.40827665) q[0];
rz(2.4093742) q[1];
sx q[1];
rz(-0.1258985) q[1];
sx q[1];
rz(0.29762038) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0746411) q[0];
sx q[0];
rz(-0.56109259) q[0];
sx q[0];
rz(-0.69343167) q[0];
x q[1];
rz(-3.1213644) q[2];
sx q[2];
rz(-2.0654404) q[2];
sx q[2];
rz(1.1927562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0875594) q[1];
sx q[1];
rz(-0.39759025) q[1];
sx q[1];
rz(1.041574) q[1];
rz(-pi) q[2];
rz(-1.4206395) q[3];
sx q[3];
rz(-0.87267002) q[3];
sx q[3];
rz(-2.0719178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97062651) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(2.4453898) q[2];
rz(-0.89037406) q[3];
sx q[3];
rz(-1.1768769) q[3];
sx q[3];
rz(-1.6740929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17906976) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(-2.9288375) q[0];
rz(0.46956024) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(2.3874217) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045005) q[0];
sx q[0];
rz(-1.8593328) q[0];
sx q[0];
rz(-2.9456784) q[0];
rz(-pi) q[1];
rz(2.5905072) q[2];
sx q[2];
rz(-2.7646825) q[2];
sx q[2];
rz(2.2152679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6282916) q[1];
sx q[1];
rz(-0.65807146) q[1];
sx q[1];
rz(0.18690303) q[1];
rz(-pi) q[2];
rz(0.92408224) q[3];
sx q[3];
rz(-1.0093186) q[3];
sx q[3];
rz(1.1738861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.135005) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(2.3593486) q[2];
rz(1.6953281) q[3];
sx q[3];
rz(-0.54636121) q[3];
sx q[3];
rz(2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049147216) q[0];
sx q[0];
rz(-0.47098422) q[0];
sx q[0];
rz(-0.96442047) q[0];
rz(1.2760705) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.6395578) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3827189) q[0];
sx q[0];
rz(-1.8028462) q[0];
sx q[0];
rz(-0.20573549) q[0];
rz(2.658074) q[2];
sx q[2];
rz(-2.3502825) q[2];
sx q[2];
rz(2.647959) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.0037476841) q[1];
sx q[1];
rz(-1.5478351) q[1];
sx q[1];
rz(2.0883043) q[1];
rz(-1.1427059) q[3];
sx q[3];
rz(-2.9799298) q[3];
sx q[3];
rz(-0.77199329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8883349) q[2];
sx q[2];
rz(-1.8586681) q[2];
sx q[2];
rz(0.9453195) q[2];
rz(0.82593289) q[3];
sx q[3];
rz(-1.632894) q[3];
sx q[3];
rz(2.6065684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(3.0655521) q[0];
sx q[0];
rz(-1.9726418) q[0];
sx q[0];
rz(-0.8031351) q[0];
rz(1.5777292) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(2.8520083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0286547) q[0];
sx q[0];
rz(-1.583033) q[0];
sx q[0];
rz(1.6763681) q[0];
x q[1];
rz(-0.8044462) q[2];
sx q[2];
rz(-0.44155332) q[2];
sx q[2];
rz(-1.785977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18145308) q[1];
sx q[1];
rz(-1.3457237) q[1];
sx q[1];
rz(-1.9378661) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2373459) q[3];
sx q[3];
rz(-1.8419918) q[3];
sx q[3];
rz(1.1734133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76558602) q[2];
sx q[2];
rz(-0.12066081) q[2];
sx q[2];
rz(0.96013367) q[2];
rz(0.58297408) q[3];
sx q[3];
rz(-2.4834902) q[3];
sx q[3];
rz(-0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4509907) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(3.1008537) q[1];
sx q[1];
rz(-2.4650885) q[1];
sx q[1];
rz(-3.0104641) q[1];
rz(-1.1004943) q[2];
sx q[2];
rz(-1.8467156) q[2];
sx q[2];
rz(2.4450532) q[2];
rz(2.1914239) q[3];
sx q[3];
rz(-0.087407268) q[3];
sx q[3];
rz(2.0711318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
