import random

# * 单工具调用 prompt 

def single_case_prompt():
    style_list = [
    '''Now finish the following task. You just need to reply with the "user query" with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the formal style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the informal style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the inspirational style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the empathetic style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the brief style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the concrete style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the humorous style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the technical style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the analytical style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the quizzical style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the pragmatic style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the assertive style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the reflective style with all parameters mentioned in it.''',
    '''Now finish the following task. You just need to reply with the "user query" in the authoritative style with all parameters mentioned in it.''',
    ]
    return \
'''Please guess the user query with the tool calling result.

Example:

Tool Information: {{
    "name": "synthesize_chemical_compound",
    "description": "Synthesizes a chemical compound based on the provided compound name.",
    "parameters": {{
        "type": "object",
        "properties": {{
            "needed_compound": {{
                "type": "string",
                "description": "The name of the chemical compound that needs to be  synthesized."
            }}
        }}
    }}
}}
Tool Calling: synthesize_chemical_compound("Asprin")
User Query: How can I synthesize Asprin?

'''+random.choice(style_list)+'''

Tool Information: {}
Tool Calling: {}
User Query: ___'''

'''The tool agent returns the tool calling according to the query user given. Please guess the user query with the tool calling now. '''

# * 多工具调用 prompt 

# * == == == == == == 

prompt_candidate_tool_selection = '''Now finish the following steps to construct an important piece of data:
1. Generate the "calling_chain": Select relevant tools to make up a proper task. The selected tools should have strong logical relationships.
2. Generate the "query": Create the query of the task based on the selected tools.
3. Generate the "reordered_calling_chain": Reorder the selected tools if needed based on the logical order of tool execution in "query".

For example:
Input:
candidate_tools = {{"daily_lib/getWeatherForecast": "Retrieve weather forecast information", "daily_lib/calculateBMI": "Calculate Body Mass Index (BMI) based on height and weight", "daily_lib/translateText": "Translate text from one language to another", "daily_lib/generateQRCode": "Generate a QR code for a given text or URL", "daily_lib/getHotelDetails": "Retrieve detailed information about a hotel", "daily_lib/getAirQualityIndex": "Retrieve the air quality index (AQI) information for a specific location", "daily_lib/searchRestaurant": "Search for a restaurant based on various criteria", "daily_lib/checkTrafficConditions": "Retrieve current traffic conditions information", "daily_lib/searchHotels": "Search for hotels based on various criteria", "daily_lib/reserveRentalCar": "Reserve a rental car for a specific location and time", "daily_lib/checkFlightAvailability": "Check the availability of flights for a specified route and date", "daily_lib/getArticleDetails": "Retrieve details of an article by providing its identifier", "daily_lib/cancelHotelReservation": "Cancel a hotel reservation", "daily_lib/callTaxi": "Request a taxi service for transportation"}}
Output:
1. calling_chain = ["daily_lib/getHotelDetails", "daily_lib/searchHotels", "daily_lib/cancelHotelReservation"]
2. query = ["Find the reserved hotel and obtain its information in order to cancel the reservation due to a schedule change."]
3. reordered_calling_chain = ["daily_lib/searchHotels", "daily_lib/getHotelDetails", "daily_lib/cancelHotelReservation"]

Please generate in the format of list [] like the example above.

Input:
1. candidate_tools = {}
Output:
1. calling_chain = [...]
2. query = [...]
3. reordered_calling_chain = [...]
'''

'''
First, generate the "calling_chain" by selecting tools with a strong logical relationship from the "candidate_tools".
Next, generate the "query" based on the calling chain.
Finally, generate the 'reordered_calling_chain' based on the tool execution order in the query.
'''

# * == == == == == == 

prompt_tool_calling_chain = '''Please generate the next tool calling according to the tool calling chain and other tool information given.
Tips:
1. You can reorder the selected tools in the tool calling chain to make it more logical.
2. To ensure the quality of the calling chain, just choose from the provided calling examples as input or make up some concrete inputs if you are proficient.
3. Try to use return values of the former tool callings as the later input parameters. Correlate individual tool callings as much as possible.

For example:
Input:
Query: Find a nearby restaurant with good reviews, then check the current traffic conditions. If the traffic is favorable, reserve a rental car for it.
Selected tools: {{"daily_lib/searchRestaurant", "daily_lib/callTaxi", "daily_lib/checkTrafficConditions"}}
Detailed tool information:
Name: daily_lib/searchRestaurant
Description: Search for a restaurant based on various criteria.
Calling examples: [["Chinese"], ["fastfood"], ["Japanese"], ["barbecue"], ["French"]]
Name: daily_lib/callTaxi
Description: Retrieve current traffic conditions information.
Calling examples: [["Central Park", "5th Avenue"], ["Union Station", "Capitol Hill"], ["Pike Place Market", "Space Needle"], ["Shibuya Crossing", "Tokyo Tower"], ["Oxford Street", "Hyde Park"]]
Name: daily_lib/checkTrafficConditions
Description: Request a taxi service for transportation.
Calling examples: [["Champs-Élysées", "morning"], ["Bourbon Street", "midnight"], ["Bolshoi Theatre", "12:00a.m."], ["Red Square", "afternoon"], ["Bolshoi Theatre", "dawn"]]
Tool calling chain:
Calling: daily_lib/searchRestaurant["Italian"]
Return: "Fortune Street, AllStars District"
Calling: daily_lib/checkTrafficConditions["Fortune Street", "afternoon"]
Return: "In good condition"
Output:
Calling: daily_lib/callTaxi["Nanjing Road", "Fortune Street"]

Now try to return the next tool calling in the brief style like the example above.

Input:
Query: {}
Selected tools: {}
Detailed tool information:
{}
Tool calling chain:
{}
Output:
Calling: ...
'''

# * == == == == == == 

# prompt_chain_to_query = '''Please generate the detailed query according to the tool information and tool calling chain.
# Tips:
# 1. Do not mention the concept of "tool" in the query.
# 2. Include all necessary parameters in the query, except for those parameters provided in the return value of previous tool callings.

# For example:
# Input:
# Tool information:
# {{"name": "daily_lib/searchRestaurant", "description": "Search for a restaurant based on various criteria", "parameters": {{"type": "object", "properties": {{"cuisine": {{"type": "str", "description": "The type of cuisine you prefer"}}}}}}}}
# {{"name": "daily_lib/checkTrafficConditions", "description": "Retrieve current traffic conditions information", "parameters": {{"type": "object", "properties": {{"location": {{"type": "str", "description": "The location for which you want to check traffic conditions"}}, "time_of_day": {{"type": "str", "description": "Specify the time of day for checking traffic conditions"}}}}}}}}
# {{"name": "daily_lib/callTaxi", "description": "Request a taxi service for transportation", "parameters": {{"type": "object", "properties": {{"pickup_location": {{"type": "str", "description": "The location where you want to be picked up"}}, "destination": {{"type": "str", "description": "The destination address where you want to go"}}}}}}}}
# Tool calling chain:
# Calling: daily_lib/searchRestaurant("Italian")
# Return: "Fortune Street, AllStars District"
# Calling: daily_lib/checkTrafficConditions("Fortune Street", "afternoon")
# Return: "In good condition"
# Calling: daily_lib/callTaxi("Nanjing Road", "Fortune Street")
# Return: {{"state": "waiting"}}
# Output:
# Query: We are on Nanjing Road and looking for a well-reviewed Italian restaurant. Please also check the current traffic conditions at the restaurant's location in the afternoon. If the traffic is favorable, proceed to reserve a rental car for us.

# Now try to generate the query like the example above.

# Input:
# Tool information:
# {}
# Tool calling chain:
# {}
# Output:
# Query: ...
# '''

def prompt_chain_to_query(): # for pymatgen
    style_list = [
    '''.''',
    ''' in the formal style.''',
    ''' in the informal style.''',
    ''' in the inspirational style.''',
    ''' in the empathetic style.''',
    ''' in the brief style.''',
    ''' in the concrete style.''',
    ''' in the humorous style.''',
    ''' in the technical style.''',
    ''' in the analytical style.''',
    ''' in the quizzical style.''',
    ''' in the pragmatic style.''',
    ''' in the assertive style.''',
    ''' in the reflective style.''',
    ''' in the authoritative style.''',
    ]
    return \
'''Please generate the detailed query according to the tool information and tool calling chain.
Tips:
1. Do not mention the concept of "tool" in the query.
2. Include all function parameters in the query. For part of them, you can express them in a different way. But all file name should be given explicitly.

For example:
Input:
Tool information:
{{"name": "daily_lib/searchRestaurant", "description": "Search for a restaurant based on various criteria", "parameters": {{"type": "object", "properties": {{"cuisine": {{"type": "str", "description": "The type of cuisine you prefer"}}}}}}}}
{{"name": "daily_lib/checkTrafficConditions", "description": "Retrieve current traffic conditions information", "parameters": {{"type": "object", "properties": {{"location": {{"type": "str", "description": "The location for which you want to check traffic conditions"}}, "time_of_day": {{"type": "str", "description": "Specify the time of day for checking traffic conditions"}}}}}}}}
{{"name": "daily_lib/callTaxi", "description": "Request a taxi service for transportation", "parameters": {{"type": "object", "properties": {{"pickup_location": {{"type": "str", "description": "The location where you want to be picked up"}}, "destination": {{"type": "str", "description": "The destination address where you want to go"}}}}}}}}
Tool calling chain:
Calling: daily_lib/searchRestaurant("Italian")
Return: "Fortune Street, AllStars District"
Calling: daily_lib/checkTrafficConditions("Fortune Street", "afternoon")
Return: "In good condition"
Calling: daily_lib/callTaxi("Nanjing Road", "Fortune Street")
Return: {{"state": "waiting"}}
Output:
Query: We are on Nanjing Road and looking for a well-reviewed Italian restaurant, please check the current traffic conditions at the restaurant's location in the afternoon. If the traffic is favorable, proceed to reserve a rental car for us.

Now try to generate the query like the example above'''+random.choice(style_list)+'''

Input:
Tool information:
{}
Tool calling chain:
{}
Output:
Query: ...
'''

# def prompt_chain_to_query():
#     style_list = [
#     '''.''',
#     ''' in the formal style.''',
#     ''' in the informal style.''',
#     ''' in the inspirational style.''',
#     ''' in the empathetic style.''',
#     ''' in the brief style.''',
#     ''' in the concrete style.''',
#     ''' in the humorous style.''',
#     ''' in the technical style.''',
#     ''' in the analytical style.''',
#     ''' in the quizzical style.''',
#     ''' in the pragmatic style.''',
#     ''' in the assertive style.''',
#     ''' in the reflective style.''',
#     ''' in the authoritative style.''',
#     ]
#     return \
# '''Please generate the detailed query according to the tool information and tool calling chain.
# Tips:
# 1. Do not mention the concept of "tool" in the query.
# 2. Include all necessary parameters in the query, except for those parameters provided in the return value of previous tool callings.

# For example:
# Input:
# Tool information:
# {{"name": "daily_lib/searchRestaurant", "description": "Search for a restaurant based on various criteria", "parameters": {{"type": "object", "properties": {{"cuisine": {{"type": "str", "description": "The type of cuisine you prefer"}}}}}}}}
# {{"name": "daily_lib/checkTrafficConditions", "description": "Retrieve current traffic conditions information", "parameters": {{"type": "object", "properties": {{"location": {{"type": "str", "description": "The location for which you want to check traffic conditions"}}, "time_of_day": {{"type": "str", "description": "Specify the time of day for checking traffic conditions"}}}}}}}}
# {{"name": "daily_lib/callTaxi", "description": "Request a taxi service for transportation", "parameters": {{"type": "object", "properties": {{"pickup_location": {{"type": "str", "description": "The location where you want to be picked up"}}, "destination": {{"type": "str", "description": "The destination address where you want to go"}}}}}}}}
# Tool calling chain:
# Calling: daily_lib/searchRestaurant("Italian")
# Return: "Fortune Street, AllStars District"
# Calling: daily_lib/checkTrafficConditions("Fortune Street", "afternoon")
# Return: "In good condition"
# Calling: daily_lib/callTaxi("Nanjing Road", "Fortune Street")
# Return: {{"state": "waiting"}}
# Output:
# Query: We are on Nanjing Road and looking for a well-reviewed Italian restaurant, please check the current traffic conditions at the restaurant's location in the afternoon. If the traffic is favorable, proceed to reserve a rental car for us.

# Now try to generate the query like the example above'''+random.choice(style_list)+'''

# Input:
# Tool information:
# {}
# Tool calling chain:
# {}
# Output:
# Query: ...
# '''

# * == == == == == == 


prompt_eval_0_shot = {
"instruction": '''Evaluate whether to invoke a tool based on the question, the existing tool-calling chain, and tool details. If a tool call is needed, return it in the exact format "Tool Calling: Tool_Name(Parameter_1, Parameter_2, ...)", with no additional text. If no tool call is needed, return only "None".

Question: 
What is the energy level of a hydrogen atom when the electron is in the second energy level?
Here is a list of tools in JSON format that you can invoke:
{{"name": "Physlib/calculateEnergyLevel", "description": "Calculates the energy of an electron in a given energy level for a specified element.", "parameters": {{"type": "object","properties": {{"element": {{"type": "str", "description": "The chemical symbol of the element."}}, "n": {{"type": "int", "description": "The principal quantum number (energy level) of the electron."}}}}}}}}

Tool Calling: Physlib/calculateEnergyLevel(element='H', n=2)

Question: 
What is the CAS number of penicillin?
Here is a list of tools in JSON format that you can invoke:
{{"name": "chemistool/getCASNumber", "description": "Retrieves the CAS number for a specified chemical compound.", "parameters": {{"type": "object", "properties": {{"compound": {{"type": "str", "description": "The common or IUPAC name of the compound."}}}}}}}}

Tool Calling: chemistool/getCASNumber(compound='penicillin')

''',
"input": '''Question:
{}{}
Here is a list of tools in JSON format that you can invoke:
{}

'''
}

# - Prompt
# prompt_eval_0_shot = {
# "instruction": '''Evaluate whether to invoke a tool based on the question, the existing tool-calling chain, and tool details. If a tool call is needed, return it in the exact format "Tool Calling: Tool_Name(Parameter_1, Parameter_2, ...)" like the example "Tool Calling: calculateEnergyLevel(element='H', n=2)", with no additional text. If no tool call is needed, return only "None".
# ''',
# "input": '''Question:
# {}
# {}
# Here is a list of tools in JSON format that you can invoke:
# {}

# Tool Calling: ...
# '''
# }

# - 第一版 Prompt eval
# prompt_eval_0_shot = {
# "instruction": '''Evaluate whether to invoke a tool based on the question, the existing tool-calling chain, and tool details. If a tool call is needed, return it in the exact format: "Tool Calling: Tool_Name(Parameter_1, Parameter_2, ...)", with no additional text. If no tool call is needed, return only "None".
# ''',
# "input": '''Question:
# {}
# {}
# Here is a list of tools in JSON format that you can invoke:
# {}

# Tool Calling: ...
# '''
# }

# 旧的single_case prompt eval_0_shot
# "Please judge whether to call a function accroding to the question, existing tool calling chain, and tool information. If you choose to return one function call in the format \"Tool Calling: Tool_Name(Parameter_1, Parameter_2, ...)\", it is imperative that no other text is included. If none function call is needed, return \"None\" only."
# '''Question:
# {}
# {}
# Here is a list of tools in JSON format that you can invoke:
# {}
# If you choose to return one function call in the format \"Tool Calling: Tool_Name(Parameter_1, Parameter_2, ...)\", it is imperative that no other text is included. If none function call is needed, return \"None\" only.
# Tool Calling: ...'''

# "Question:{}\nHere is a list of tools in JSON format that you can invoke:\n{}\nOnly need to return the function call in the format like Tool_Name(Parameter_1, Parameter_2, ...), NO other text MUST be included."

# *** 

prompt_response_generation = """Please anwser the question in one sentence based on the given tool calling result.

Question:
{}

Tool calling chain:
{}

Answer:
"""

# for GPT GOLD ANSWER
# prompt_response_generation = """Please answer the question concisely based on the tool calling chain, other information given and your knowledge.

# Question: {}

# Tool calling chain:
# {}

# Response:"""

prompt_response_evaluation = """Please evaluate whether the given answer to the question is correct according to the gold answer. Return "Yes" if it is correct, otherwise return "No".

Question:
{}

Gold Answer:
{}

Given Answer:
{}

Whether the given answer is correct:
"""

prompt_answer = '''Please answer the question in one sentence.

Question:
{}

Answer:
'''