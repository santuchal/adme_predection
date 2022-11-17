## ADMET prediction using ML

This post is a step-by-step implementation of my approach to calculating the ADME(T) properties and it is meant for sharing, studying, and critiquing by fellow researchers who are new and interested in this topic. Most of the research paper that is implemented can be found in **ref** directory. 

### **Physicochemical Properties:**

*   [x] **Molecular Weight**
*   [ ] **Volume**
*   [ ] **Density**
*   [x] **Heavy Atoms**
*   [x] **Aromatic Heavy Atoms**
*   [x] **Fraction Csp3**
*   [x] **Rotatable Bonds**
*   [x] **Hetro Atoms**
*   [x] **H-Bond Acceptors**
*   [x] **H-Bond Donors**
*   [x] **Ring Count**
*   [x] **Aromatic Ring Count**
*   [x] **Stereo Centers**
*   [x] **Molar Refractivity**
*   [x] **tPSA**
*   [x] **LogP**
*   [ ] **pKa**
*   [x] **LogS**
*   [ ] **LogD7.4**

---

***LogS*** 



The work of Delaney provided the following linear regression equation:

    LogS = 0.16 - 0.63 cLogP - 0.0062 MW + 0.066 RB - 0.74 AP

The reproduction by Pat Walters

provided the following:

    LogS = 0.26 - 0.74 LogP - 0.0066 MW + 0.0034 RB - 0.42 AP


---

 #ADME #adme #Physicochemical